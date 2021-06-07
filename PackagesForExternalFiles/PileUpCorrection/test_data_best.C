///////////////////////////////////////////////////////////////////////////
//  This macro performs unfolding to extract true refmult3 distribution 
//  from sinle-collision events, based on experimental measured distribution 
//  including pilupe events.
//
//  usege:
//  root -l test_data_best.C++
//
//  for error estimation:
//  root -l test_data_best.C++'(1,iBS)'
//  , where iBS is integer number from 0-99.
//
//  One can repeat the unfolding by varying pileup probability, alpha, to 
//  determine the best parameter which can describe the data.
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TProfile.h>

using std::cout;
using std::endl;

const double alpha = 0.00556325; // best
//const double alpha = 0.00467747; // best - sigma
//const double alpha = 0.00649443; // best + sigma
const int nitr = 100;
const double sc = 0.5;

const int nbin = 200;
const int FitStart = 0;
const int FitEnd = 150;

void Diff( TH1* hData, TH1 *hMC ){
	TH1 *hMCtmp = (TH1D*)hMC->Clone("hMCtmp");
	hMC->Reset();
	for(int ibin=FitStart; ibin<FitEnd; ibin++){
		const double cont1 = hData->GetBinContent(ibin+1);
		const double cont2 = hMCtmp->GetBinContent(ibin+1);
		hMC->SetBinContent(ibin+1,cont1-cont2);
	}
}


void test_data_best( const bool IsBootstrap=false, const int iBS=0 ){

	gRandom->SetSeed(0);

	cout <<"Working on alpha="<< alpha << endl;

	// Experimental data 
	TFile *fData = TFile::Open("DataFromSam.root","READ");
	TH1 *hData_all = (TH1D*)fData->Get("FXTMult_noProtons_orAnti");
	TH1 *hData_all_bs = (TH1D*)hData_all->Clone("hData_all_bs");
	const double nevent_data = hData_all->GetEntries();
	if(IsBootstrap){
		hData_all->Reset();
		for(int ievent=0; ievent<nevent_data; ievent++){
			hData_all->Fill((int)(hData_all_bs->GetRandom()+0.5));
		}
	}
	hData_all->Scale(nevent_data);
	cout << nevent_data <<" events found in data."<< endl;

	// Initial condition of MC distribution.
	// This can be any probability distribution like uniform distribution
	TFile *fMC = TFile::Open("glauber_pu1per_eff.root","READ");
	TH1 *hMC_true_tmp = (TH1D*)fMC->Get("hMult_true");
	TH1 *hMC_all_tmp = (TH1D*)fMC->Get("hMult");
	TH1 *hMC_true = new TH1D("hMC_true","hMC_true",200,-0.5,199.5);
	TH1 *hMC_all = new TH1D("hMC_all","hMC_all",200,-0.5,199.5);
	for(int ibin=0; ibin<200; ibin++){
		hMC_true->Fill(ibin,hMC_true_tmp->GetBinContent(ibin+1));
		hMC_all->Fill(ibin,hMC_all_tmp->GetBinContent(ibin+1));
	}
	const double nevent_mc = hMC_all->Integral();
	cout << nevent_mc <<" events found in MC."<< endl;
	const double scale_data = hData_all->Integral(FitStart,FitEnd);
	const double scale_mc = hMC_all->Integral(FitStart,FitEnd);

	hData_all->Sumw2();
	hData_all->Scale(1/scale_data);
	hMC_true->Sumw2();
	hMC_all->Sumw2();
	hMC_true->Scale(1/scale_mc);
	hMC_all->Scale(1/scale_mc);


	int N_all, Ntrue, Npu;
	TH1 *hMC_true_itr[nitr+1];
	TH1 *hMC_all_itr[nitr+1];
	TH1 *hRM[nbin][nitr+1];
	for(int iitr=0; iitr<nitr+1; iitr++){
		hMC_true_itr[iitr] = (TH1D*)hMC_true->Clone(Form("hMC_true_itr_%d",iitr));
		hMC_all_itr[iitr] = (TH1D*)hMC_all->Clone(Form("hMC_all_itr_%d",iitr));
		hMC_true_itr[iitr]->Reset();
		hMC_all_itr[iitr]->Reset();
		for(int ibin=0; ibin<nbin; ibin++){
			hRM[ibin][iitr] = (TH1D*)hMC_all->Clone(Form("hRM_%d_%d",ibin,iitr));
			hRM[ibin][iitr]->Reset();
		}
	}
	TH1 *hMC_all_diff_itr[nitr+1];
	TH1 *hMC_true_diff_itr[nitr+1];
	TH1 *hMC_tmp[nitr+1];
	TH2 *hRM_final = new TH2D("hRM_final","hRM_final",nbin,-0.5,nbin-0.5,nbin,-0.5,nbin-0.5);
	TH1 *hMC_all_final = new TH1D("hMC_all_final","hMC_all_final",nbin,-0.5,nbin-0.5);
	// Iteration
	for(int iitr=0; iitr<nitr; iitr++){
		if(iitr==0){
			hMC_tmp[iitr] = (TH1D*)hMC_true->Clone(Form("hMC_tmp_%d",iitr));
		}
		else {
			hMC_tmp[iitr] = (TH1D*)hMC_true_itr[iitr]->Clone(Form("hMC_tmp_%d",iitr));
		}
		hMC_true_itr[iitr]->Reset();
		for(int ievent=0; ievent<nevent_mc; ievent++){
			Ntrue = (int)(hMC_tmp[iitr]->GetRandom()+0.5);
			N_all = Ntrue;
			if(gRandom->Rndm()<alpha){
				Npu = (int)(hMC_tmp[iitr]->GetRandom()+0.5);
				N_all += Npu; 
				if(N_all>=nbin) continue;
				hRM[N_all][iitr]->Fill(Ntrue);
				hRM[N_all][iitr]->Fill(Npu);
				if(iitr==nitr-1) hRM_final->Fill(Ntrue,Npu);
			}
			else{
				hMC_true_itr[iitr]->Fill(N_all);
			}
			hMC_all_itr[iitr]->Fill(N_all);
			if(iitr==nitr-1) hMC_all_final->Fill(N_all);
		}
		hMC_true_itr[iitr]->Scale(1/hMC_all_itr[iitr]->Integral(FitStart,FitEnd));
		hMC_all_itr[iitr]->Scale(1/hMC_all_itr[iitr]->Integral(FitStart,FitEnd));
		hMC_all_diff_itr[iitr] = (TH1D*)hMC_all_itr[iitr]->Clone(Form("hMC_all_diff_itr_%d",iitr));
		Diff(hData_all,hMC_all_diff_itr[iitr]);
		for(int ibin=0; ibin<nbin; ibin++){
			const double nevent_bin = hRM[ibin][iitr]->Integral();
			if(nevent_bin==0) continue; 
			hRM[ibin][iitr]->Scale(1/nevent_bin);
		}
		//// Apply response matrix to convert to the true coordinate
		hMC_true_diff_itr[iitr] = (TH1D*)hMC_all_diff_itr[iitr]->Clone(Form("hMC_true_diff_itr_%d",iitr));
		hMC_true_diff_itr[iitr]->Reset();
		for(int ibin=0; ibin<nbin; ibin++){
			const double diff = hMC_all_diff_itr[iitr]->GetBinContent(ibin+1);
			if(diff==0) continue;
			for(int jbin=0; jbin<nbin; jbin++){
				if(hRM[ibin][iitr]->Integral()==0) continue;
				const double cont = hRM[ibin][iitr]->GetBinContent(jbin+1);
				if(cont==0) continue;
				hMC_true_diff_itr[iitr]->Fill(jbin,cont*diff);
			//	cout << ibin <<"  "<< cont*diff << endl;
			}
		}
		for(int ibin=0; ibin<nbin; ibin++){
			hMC_true_diff_itr[iitr]->SetBinError(ibin+1,0);
			hMC_true_itr[iitr]->SetBinError(ibin+1,0);
			hMC_all_itr[iitr]->SetBinError(ibin+1,0);
		}
		//// Apply correction functions
		for(int ibin=0; ibin<nbin; ibin++){
			const double init = hMC_true_itr[iitr]->GetBinContent(ibin+1);
			const double cor = sc*hMC_true_diff_itr[iitr]->GetBinContent(ibin+1);
			if(init+cor>0){
				hMC_true_itr[iitr+1]->SetBinContent(ibin+1,init+cor);
			//	cout << ibin <<"  "<< init+cor << endl;
			}
		}
		cout <<"Working on iteration #"<< iitr << endl;
	}


	TCanvas *c1 = new TCanvas("c1","c1",0,0,1400,600);
	c1->Divide(2,1);
	c1->cd(1);
	gPad->SetLogy();
	hData_all->Draw();
	hData_all->GetXaxis()->SetRangeUser(0,170);
	hMC_all_itr[0]->Draw("same");
	hMC_all_itr[1]->Draw("same");
	hMC_all_itr[nitr-1]->Draw("same");
	hMC_all_itr[0]->SetMarkerColor(kRed);
	hMC_all_itr[0]->SetLineColor(kRed);
	hMC_all_itr[1]->SetMarkerColor(kBlue);
	hMC_all_itr[1]->SetLineColor(kBlue);
	hMC_all_itr[nitr-1]->SetMarkerColor(kGreen+1);
	hMC_all_itr[nitr-1]->SetLineColor(kGreen+1);
	c1->cd(2);
	gPad->SetLogy();
	hMC_true_itr[0]->Draw("same");
	hMC_true_itr[1]->Draw("same");
	hMC_true_itr[0]->GetXaxis()->SetRangeUser(0,170);
	hMC_true_itr[nitr-1]->Draw("same");
	hMC_true_itr[0]->SetMarkerColor(kRed);
	hMC_true_itr[0]->SetLineColor(kRed);
	hMC_true_itr[1]->SetMarkerColor(kBlue);
	hMC_true_itr[1]->SetLineColor(kBlue);
	hMC_true_itr[nitr-1]->SetMarkerColor(kGreen+1);
	hMC_true_itr[nitr-1]->SetLineColor(kGreen+1);


	TCanvas *c2 = new TCanvas("c2","c2",0,0,1200,900);
	c2->Divide(2,2);
	c2->cd(3);
	gPad->SetLogy();
	hRM[50][0]->Draw();
	c2->cd(1);
	hMC_all_diff_itr[0]->Draw();
	hMC_all_diff_itr[1]->Draw("same");
	hMC_all_diff_itr[2]->Draw("same");
	hMC_all_diff_itr[3]->Draw("same");
	hMC_all_diff_itr[nitr-1]->Draw("same");
	hMC_all_diff_itr[0]->GetXaxis()->SetRangeUser(0,170);
	hMC_all_diff_itr[0]->SetMinimum(-0.005);
	hMC_all_diff_itr[0]->SetMaximum(0.005);
	hMC_all_diff_itr[0]->SetMarkerStyle(24);
	hMC_all_diff_itr[0]->SetMarkerColor(kRed);
	hMC_all_diff_itr[0]->SetLineColor(kRed);
	hMC_all_diff_itr[1]->SetMarkerStyle(24);
	hMC_all_diff_itr[1]->SetMarkerColor(kBlue);
	hMC_all_diff_itr[1]->SetLineColor(kBlue);
	hMC_all_diff_itr[2]->SetMarkerStyle(24);
	hMC_all_diff_itr[2]->SetMarkerColor(kGreen+2);
	hMC_all_diff_itr[2]->SetLineColor(kGreen+2);
	hMC_all_diff_itr[3]->SetMarkerStyle(24);
	hMC_all_diff_itr[3]->SetMarkerColor(kMagenta);
	hMC_all_diff_itr[3]->SetLineColor(kMagenta);
	hMC_all_diff_itr[nitr-1]->SetMarkerStyle(24);
	hMC_all_diff_itr[nitr-1]->SetMarkerColor(kBlack);
	hMC_all_diff_itr[nitr-1]->SetLineColor(kBlack);

	c2->cd(2);
	hMC_true_diff_itr[0]->Draw();
	hMC_true_diff_itr[1]->Draw("same");
	hMC_true_diff_itr[2]->Draw("same");
	hMC_true_diff_itr[3]->Draw("same");
	hMC_true_diff_itr[nitr-1]->Draw("same");
	hMC_true_diff_itr[0]->GetXaxis()->SetRangeUser(0,170);
	hMC_true_diff_itr[0]->SetMinimum(-0.005);
	hMC_true_diff_itr[0]->SetMaximum(0.005);
	hMC_true_diff_itr[0]->SetMarkerStyle(24);
	hMC_true_diff_itr[0]->SetMarkerColor(kRed);
	hMC_true_diff_itr[0]->SetLineColor(kRed);
	hMC_true_diff_itr[1]->SetMarkerStyle(24);
	hMC_true_diff_itr[1]->SetMarkerColor(kBlue);
	hMC_true_diff_itr[1]->SetLineColor(kBlue);
	hMC_true_diff_itr[2]->SetMarkerStyle(24);
	hMC_true_diff_itr[2]->SetMarkerColor(kGreen+2);
	hMC_true_diff_itr[2]->SetLineColor(kGreen+2);
	hMC_true_diff_itr[3]->SetMarkerStyle(24);
	hMC_true_diff_itr[3]->SetMarkerColor(kMagenta);
	hMC_true_diff_itr[3]->SetLineColor(kMagenta);
	hMC_true_diff_itr[nitr-1]->SetMarkerStyle(24);
	hMC_true_diff_itr[nitr-1]->SetMarkerColor(kBlack);
	hMC_true_diff_itr[nitr-1]->SetLineColor(kBlack);

	TProfile *pMC_true_itr = new TProfile("pMC_true_itr","pMC_true_itr",nbin,0,nbin,"S");
	TProfile *pMC_all_itr = new TProfile("pMC_all_itr","pMC_all_itr",nbin,0,nbin,"S");
	for(int ibin=0; ibin<nbin; ibin++){
		pMC_true_itr->Fill(ibin,hMC_true_itr[nitr-1]->GetBinContent(ibin+1));
		pMC_all_itr->Fill(ibin,hMC_all_itr[nitr-1]->GetBinContent(ibin+1));
	}

	TFile *fOut;
	if(IsBootstrap)  fOut = new TFile(Form("test_data_best_bs%d.root",iBS),"RECREATE");
	if(!IsBootstrap) fOut = new TFile("test_data_best.root","RECREATE");
//	if(IsBootstrap)  fOut = new TFile(Form("test_data_best_minus_bs%d.root",iBS),"RECREATE");
//	if(!IsBootstrap) fOut = new TFile("test_data_best_minus.root","RECREATE");
//	if(IsBootstrap)  fOut = new TFile(Form("test_data_best_plus_bs%d.root",iBS),"RECREATE");
//	if(!IsBootstrap) fOut = new TFile("test_data_best_plus.root","RECREATE");
	for(int iitr=0; iitr<nitr; iitr++){
		hMC_true_itr[iitr]->Scale(hMC_true_itr[iitr]->GetEntries());
		hMC_all_itr[iitr]->Scale(hMC_all_itr[iitr]->GetEntries());
		hMC_true_itr[iitr]->Write();
		hMC_all_itr[iitr]->Write();
		hMC_true_diff_itr[iitr]->Write();
		hMC_all_diff_itr[iitr]->Write();
	}
	pMC_true_itr->Write();
	pMC_all_itr->Write();
	hData_all->SetName("hData_all");
	hData_all->Write();
	hMC_true->SetName("hMC_true");
	hMC_all->SetName("hMC_all");
	hMC_true->Write();
	hMC_all->Write();
	hMC_all_final->Write();
	hRM_final->Write();

	fOut->Write();
	fOut->Close();

}
