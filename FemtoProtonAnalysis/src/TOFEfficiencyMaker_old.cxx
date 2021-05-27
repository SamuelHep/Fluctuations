#include <iostream>
#include <vector>
#include <utility>
#include <numeric>

#include "TMath.h"
#include "TFile.h"
#include "TH3D.h"
#include "TChain.h"
#include "TRandom3.h"

#include "StFemtoTrack.h"
#include "StFemtoEvent.h"

#include "TGraphAsymmErrors.h"
#include "InputParameterList.h"
#include "CumulantProfileContainer.h"
#include "MomentFunctions.h"
#include "ProtonEfficiency.h"

#include "analysisUtil.h"

#include "TOFEfficiencyMaker.h"

TOFEfficiencyMaker::TOFEfficiencyMaker() : MAXFXT3(200), PTBINS(200), ETABINS(460), mass(0.938272)
{

 event = new StFemtoEvent();  

  
 std::vector<TString> h_names = {"match_nhf_10","match_nhf_12","match_nhf_15","match_m2_plus0.05","match_m2_minus0.05","all_nhf_10","all_nhf_12","all_nhf_15","all","all_pos","all_NSigCut","all_dca"};

 hist_names = h_names; 

}

void TOFEfficiencyMaker::LoopEvents(TChain * tc,
				    long int entries,
				    InputParameterList & pl,
				    TString outfilename
				    )
{
 
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  for ( auto &name : hist_names ) TOFHists[ name ] = new TH3D(name,"", MAXFXT3, -0.5,MAXFXT3-0.5, PTBINS, 0.2, 2.2, ETABINS, 0, 2.3); 
  
  for (int iEntry=0;iEntry<entries;iEntry++)
    {
      
      //      cout << iEntry << endl;
      tc->GetEntry(iEntry);
      
      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      
      //Make Event Cuts
      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;
           
      //Begin Track Loop
      for (int iTrack=0;iTrack<event->GetEntries();iTrack++)
	{
	  StFemtoTrack trk = event->GetFemtoTrack(iTrack);
	  
	  //Get Track Variables
	  double pt = trk.GetPt();
	  double pz = trk.GetPz();
	  short nhitsfit = trk.GetNHitsFit();
	  short charge = trk.GetCharge();
	  double dcaX  = trk.GetDcaX();
	  double dcaY  = trk.GetDcaY();
	  double dcaZ  = trk.GetDcaZ();
	  double dca = sqrt( dcaX*dcaX + dcaY*dcaY + dcaZ*dcaZ );

	  TVector3 pMom = TVector3(trk.GetPx(),trk.GetPy(),trk.GetPz());
	  double p = pMom.Mag();

	  double nsigpro = trk.GetNSigmaProton() - getProtonDedxMeanShift(p);
	  double nsigpi = trk.GetNSigmaPion();
	  double tofmass = trk.GetTofMass();
	  double tofmass2 = (tofmass < 0 ) ?  -1 : tofmass*tofmass;
	  
	  //calculate the energy and rapidity
	  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(mass,2) );
	  Double_t rap = 0.5 * log( ( energy + pz) / (energy - pz) );
	  rap = fabs(rap);
	  Double_t eta = 0.5*TMath::Log( ( p + pz )/( p - pz ) );
	  eta = fabs(eta);

	  //Make Track Cuts
	  TOFHists[ "all" ]->Fill( fxt3, pt, eta ); 
	  if (charge<0)   continue;
	  TOFHists[ "all_pos" ]->Fill( fxt3, pt, eta ); 
	  if ( fabs(nsigpro) > pl.NSigmaProtonCut()) continue;
	  TOFHists[ "all_NSigCut" ]->Fill( fxt3, pt, eta ); 
	  if ( dca > pl.Dca() ) continue;
	  TOFHists[ "all_dca" ]->Fill( fxt3, pt, eta ); 


	  bool nhitsfitGood[3] = {false,false,false};
	  if ( nhitsfit >= 10) nhitsfitGood[0] = true;  
	  if ( nhitsfit >= 12) nhitsfitGood[1] = true;  
	  if ( nhitsfit >= 15) nhitsfitGood[2] = true;  

	  //if pass all cuts, make a pair of the charge and the effficiency
	  double trackEff = -1;
	  bool tofmatch[3] ={false,false,false};
	  if ( tofmass2 > 0.6 && tofmass2 < 1.2) tofmatch[0] =true;
	  if ( tofmass2 > 0.65 && tofmass2 < 1.25) tofmatch[1] =true;
	  if ( tofmass2 > 0.55 && tofmass2 < 1.15) tofmatch[2] =true;

	  //All TPC Tracks
	  if ( nhitsfitGood[0] ) TOFHists[ "all_nhf_10" ]->Fill( fxt3, pt, eta ); 
	  if ( nhitsfitGood[1] ) TOFHists[ "all_nhf_12" ]->Fill( fxt3, pt, eta ); 
	  if ( nhitsfitGood[2] ) TOFHists[ "all_nhf_15" ]->Fill( fxt3, pt, eta ); 

	  //TOF Match histograms
	  if ( tofmatch[0] )
	    {
	      if ( nhitsfitGood[0] ) TOFHists[ "match_nhf_10" ]->Fill( fxt3, pt, eta ); 
	      if ( nhitsfitGood[1] ) TOFHists[ "match_nhf_12" ]->Fill( fxt3, pt, eta ); 
	      if ( nhitsfitGood[2] ) TOFHists[ "match_nhf_15" ]->Fill( fxt3, pt, eta ); 
	    }

	  if ( tofmatch[1] && nhitsfitGood[0] ) TOFHists[ "match_m2_plus0.05" ]->Fill( fxt3, pt, eta ); 
	  if ( tofmatch[2] && nhitsfitGood[0] ) TOFHists[ "match_m2_minus0.05" ]->Fill( fxt3, pt, eta ); 

	}
    }


  TFile * outfile = new TFile( outfilename, "recreate");
  outfile->cd();

  for ( auto &hist_map : TOFHists )
    {
      hist_map.second->Write();
    }

  outfile->Close();

}


void TOFEfficiencyMaker::StudyTOFEfficiency(TString filename,TString outfilename)
{


  //CURRENT CENTRALITY DEFINITION
  std::vector<double> binLabels = { 47, 70, 107, 157, 219, 282,  326};
  std::vector<int> binEdges     = { 4,    6,    10,    16,   25,  38, 48, 79};


  auto file = new TFile(filename,"read");
  for ( auto &name : hist_names ) TOFHists[ name ] = (TH3D*) file->Get(name);
  
  std::vector<TH2D*> ratioHists;
  std::vector<std::pair < TString , TString >> match_all_pair;
  match_all_pair.push_back( std::make_pair("match_nhf_10","all_nhf_10"));
  //  match_all_pair.push_back( std::make_pair("match_nhf_12","all_nhf_12"));
  //  match_all_pair.push_back( std::make_pair("match_nhf_15","all_nhf_15"));
  //  match_all_pair.push_back( std::make_pair("match_m2_plus0.05","all_nhf_10"));
  //  match_all_pair.push_back( std::make_pair("match_m2_minus0.05","all_nhf_10"));
  
  TString par = "pol(12)";
  std::vector<TCanvas*> can_vec;
  std::vector<TF1*> function_vec;
  std::vector<TGraph*> graph_vec;

  for ( auto &pair_iter : match_all_pair )
    {
      TString ratioName = pair_iter.first;
      ratioName.ReplaceAll("match_","");

      TH3D * h3_match = TOFHists[ pair_iter.first ];
      TH3D * h3_all = TOFHists[ pair_iter.second ];
      
      int minBin = h3_all->GetXaxis()->FindBin( binEdges[0] );
      int maxBin = h3_all->GetXaxis()->FindBin( binEdges[7] );
      
      h3_match->GetXaxis()->SetRange( minBin, maxBin);
      h3_all->GetXaxis()->SetRange( minBin, maxBin);
      TH2D * h2_match = (TH2D*) h3_match->Project3D("yz");
      TH2D * h2_all   = (TH2D*) h3_all->Project3D("yz");
    
      int h2_all_bins = h2_all->GetNbinsX();
      
      std::vector< TGraphAsymmErrors* > gr_vec(h2_all_bins,(TGraphAsymmErrors*) NULL);
      std::vector< TF1* > f_vec(h2_all_bins,(TF1*) NULL);

      for (int i=0;i<gr.size();i++)
	{
	  TH1D * pass = h2_match->ProjectionY(TString::Format("pass_%s_%i",ratioName.Data(),i ),i,i );
	  TH1D * total = h2_all->ProjectionY(TString::Format("total_%s_%i",ratioName.Data(),i ),i,i );

	  if ( pass->GetEntries() == 0 ) continue;
	  
	  gr_vec[i] = new TGraphAssymErrors(pass,total,"cp");
	  gr_vec[i]->SetName(TString::Format("ratioGraph_%s_%i",ratioName.Data(),i ));
	  f_vec[i] = =new TF1(TString::Format("ratioFunc_%s_%i",ratioName.Data(),i ),par,0.2,2.2);	  
	}

      std::vector<TCanvas*> canvases = PlotSlices( f_vec, gr_vec,ratioName);
      
      can_vec.insert( can_vec.end(), canvases.begin(),canvases.end());
      graph_vec.insert( graph_vec.end(), gr_vec.begin(), gr_vec.end());
      function_vec.insert( function_vec.end(), f_vec.begin(), f_vec.end());

    }
      
    
  TFile * outfile = new TFile(outfilename,"recreate");
  outfile->cd();
  outfile->mkdir("Canvases");
  outfile->cd("Canvases");
  for ( auto &h : can_vec )
    {
      h->Write();
    }

  outfile->cd("..");
  outfile->mkdir("Graphs");
  outfile->cd("Graphs");
  for ( auto &h : graph_vec )
    {
      h->Write();
    }

  outfile->cd("..");
  outfile->mkdir("Fits");
  outfile->cd("Fits");
  for ( auto &h : function_vec )
    {
      h->Write();
    }


}

std::vector<TCanvas*> TOFEfficiencyMaker::PlotSlices(std::vector<TF1*> &f, std::vector<TGraphAsymmErrors> &gr, TString canName)
{

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  int nSlices =0; 
  for ( int i=0;i<f.size();i++)
    {
      if ( !f || !gr ) continue;
      nSlices++;
    }

  int nPerPage = 9;
  int nPages = nSlices/nPerPage;
  if ( nPages % nPerPage != 0 ) nPages++;
  
  cout << "nSlices  = " << nSlices << endl;
  cout << "nPerPage = " << nPerPage << endl;
  cout << "nPages   = " << nPages << endl;

  std::vector<TCanvas*> can_vec(nPages,(TCanvas*) NULL);
  int canIndex=0;

  for ( auto &v : can_vec)
    {
      v = new TCanvas(TString::Format("can_%s_%i",canName.Data(),canIndex++),"",1200,900);
      v->Divide(3,3,0.0,0.0);
    }

  canIndex=0;
  int padIndex=1;

  for ( int i=2;i<f1_vec.size();i++)
    {

      TF1  * fit  = f[i];
      TGraphAsymmErrors * graph = gr[i];
      if ( !graph || !fit ) continue;

      can_vec[ canIndex ]->cd( padIndex );
      if ( padIndex == 3 || padIndex == 6 || padIndex == 9 ) gPad->SetRightMargin(0.02);
      if ( padIndex > 6 )
	{ 
	  graph->GetYaxis()->SetRangeUser(0,1.19);
	  graph->GetXaxis()->SetRangeUser(0.001,2.45);
	}
      else
	{
	  graph->GetYaxis()->SetRangeUser(0.001,1.19);
	  graph->GetXaxis()->SetRangeUser(0.001,2.45);
	}

      //      gPad->SetGrid(1,1);
      gPad->SetTicks(1,1);
      graph->SetMarkerStyle( kFullCircle );
      graph->SetMarkerSize( 0.5 );
      graph->GetXaxis()->SetTitle("p_{T} (GeV)");
      graph->GetYaxis()->SetTitle("Eff.");
      //      graph->SetStats(false);
      graph->Draw("ap");
      fit->SetLineStyle(9);
      fit->SetLineColor(kRed);
      fit->Draw("same");

      TString pars1 = TString::Format("Fit Parameters: [0]=%.3f, [1]=%.3f, [2]=%.3f",fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2));
      //      TString pars2 = TString::Format("[3]=%02f, [4]=%02f, [5]=%02f",fit->GetParameter(3),fit->GetParameter(4),fit->GetParameter(5));
      TString pars2 = TString::Format("                [3]=%.3f, [4]=%.3f",fit->GetParameter(3),fit->GetParameter(4));

      TLatex * text = new TLatex();
      text->SetTextSize(0.046);
      //      text->DrawLatexNDC(0.2,0.25, pars1 );
      //      text->DrawLatexNDC(0.2,0.15, pars2 );
      text->DrawLatexNDC(0.2,0.35, fit_name[ name_vec[i] ] );


      padIndex++;
      if ( padIndex == (nPerPage+1) )
	{
	  padIndex=1;
	  canIndex++;
	}     
    }

  return can_vec;

}
