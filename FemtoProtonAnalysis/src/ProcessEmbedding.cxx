#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "ProcessEmbedding.h"


using namespace std;

ProcessEmbedding::ProcessEmbedding(TString embed_filename,TString outfilename, TString particleName, InputParameterList & pl,TString sysName)
{

  gStyle->SetOptStat(0);


  _sysName = sysName;
  _pl = pl;

  if ( particleName = "p" ) particle_mass = 0.938272;
  if ( particleName = "pi" ) particle_mass = 0.13957;
  if ( particleName = "K" ) particle_mass = 0.497672;

  std::vector<int> bin_edges     = { 4, 6, 10, 16, 25,  38, 48, 79};

  _bin_edges  = bin_edges;

  TFile * f = new TFile(embed_filename,"read");
  tree = (TTree*) f->Get("EmbeddingTracks");

  tree->SetBranchAddress( "FXTMult3" , &FXTMult3 ); 
  tree->SetBranchAddress( "vtx" , &vtx ); 
  tree->SetBranchAddress( "vty" , &vty ); 
  tree->SetBranchAddress( "vtz" , &vtz ); 
  tree->SetBranchAddress( "rcPt" , &rcPt ); 
  tree->SetBranchAddress( "rcP" , &rcP ); 
  tree->SetBranchAddress( "rcEta" , &rcEta ); 
  tree->SetBranchAddress( "rcPhi" , &rcPhi ); 
  tree->SetBranchAddress( "rcCharge" , &rcCharge ); 
  tree->SetBranchAddress( "rcNfit" , &rcNfit ); 
  tree->SetBranchAddress( "rcNposs" , &rcNposs ); 
  tree->SetBranchAddress( "rcDca" , &rcDca ); 
  tree->SetBranchAddress( "rcDcaXy" , &rcDcaXy ); 
  tree->SetBranchAddress( "rcDcaZ" , &rcDcaZ );
  tree->SetBranchAddress( "mcPt" , &mcPt ); 
  tree->SetBranchAddress( "mcP" , &mcP ); 
  tree->SetBranchAddress( "mcEta" , &mcEta ); 
  
  InitHistograms();
  Run();
  RatioHistograms();
  FitAllCentSlices();
  cout << "Finshed Fit" << endl;
  PlotSlices();
  
  //Write everything to file
  TFile * outFile = new TFile( outfilename,"recreate");

  h_ratio->Write();
  h_emb->Write();
  h_reco->Write();

  for ( auto &v : h1 ) v.second->Write();
  for ( auto &v : f1 ) v.second->Write();

  outFile->cd("..");
  outFile->mkdir("Canvases");
  outFile->cd("Canvases");
  for ( auto &can : _canvases ) can->Write();

  outFile->Close();

  for ( auto &v : h1 ) delete v.second;
  for ( auto &v : f1 ) delete v.second;
  for ( auto &can : _canvases ) delete can;

}

ProcessEmbedding::~ProcessEmbedding()
{
  
  

}

Int_t ProcessEmbedding::Run()
{

  double proton_mass = 0.938272; // GeV

  for (long int i=0;i<tree->GetEntries(); i++)
    {
      tree->GetEntry(i);

      double rapidity = fabs(Rapidity( mcEta, mcPt, particle_mass ));

      if ( FXTMult3 < 4 ) continue;
      if ( FXTMult3 > 79 ) continue;
      //      if ( rcCharge < 0 ) continue;
			    
      h_emb->Fill( rapidity, mcPt );
      if ( rcNfit < _pl.NHitsFitMin() ) continue;
      if ( (double)rcNfit/(double)rcNposs < 0.51 ) continue;
      if ( rcDca > _pl.Dca() ) continue;

      if ( rcPt > 0.0 ) h_reco->Fill( rapidity, mcPt );
    }

}


Int_t ProcessEmbedding::InitHistograms()
{
  
  h_emb = new TH2D(  "emb", "", 110, 0, 2.2, 150, 0.2, 3.2 ); 
  h_reco = new TH2D(  "reco", "", 110, 0, 2.2, 150, 0.2, 3.2 ); 

  return 0;
}

TGraphAsymmErrors * ProcessEmbedding::RatioSlice(TString name, int ibin)
{
  TH1D * h1_reco = h_reco->ProjectionY(TString::Format("reco_%s",name.Data()),ibin,ibin);
  TH1D * h1_emb = h_emb->ProjectionY(TString::Format("emb_%s",name.Data()),ibin,ibin);
  TGraphAsymmErrors * tga = new TGraphAsymmErrors( h1_reco, h1_emb, "cp" );
  tga->SetName(TString::Format("Ratio_Sys_%s_%s",_sysName.Data(),name.Data()));
  return tga;
}

Int_t ProcessEmbedding::RatioHistograms()
{

  h_ratio = (TH2D*) h_reco->Clone();
  h_ratio->SetName( "ratio" );
  h_ratio->Sumw2();
  h_emb->Sumw2();
  h_ratio->Divide( h_emb );

  return 0;
}


Double_t ProcessEmbedding::Rapidity(double eta, double pT, double mass)
{
  return TMath::Log( ( sqrt( pow(mass,2.) + pow(pT,2.)*pow(TMath::CosH( eta ),2.) ) + pT*TMath::SinH( eta ) ) / sqrt( pow(mass,2.) + pow(pT,2.) ) );
}

void ProcessEmbedding::FitAllCentSlices()
{

  for ( int ibin=0; ibin< h_ratio->GetNbinsX(); ibin++)
    {
      TString name = TString::Format("slice_bin%i",ibin);
      TString yName = TString::Format("Rapidity range [%.2f,%.2f]",h_ratio->GetXaxis()->GetBinLowEdge(ibin),h_ratio->GetXaxis()->GetBinLowEdge(ibin)+h_ratio->GetXaxis()->GetBinWidth(ibin));
      
      h1[ name ] = RatioSlice( name, ibin );  //(TH1D*) h_ratio->ProjectionY( name , ibin, ibin );
      if ( h1[ name ]->GetN() == 0 ) continue;
      f1[ name ] = new TF1( TString::Format("f_Sys_%s_%s",_sysName.Data(),name.Data()), "[0]*exp(-1.0*[1]/(x^[2]))+[3]*x+[4]*x*x",0.1,3.2);

      f1[ name ]->SetParameter(0,0.95);
      f1[ name ]->SetParLimits(0,0.1,1.0);

      f1[ name ]->SetParameter(1,0.2);
      f1[ name ]->SetParLimits(1,0.00001,1.0);

      f1[ name ]->SetParameter(2,2.0);
      f1[ name ]->SetParLimits(2,0.1,50.0);

      f1[ name ]->SetParameter(3,0.0);	
      f1[ name ]->SetParameter(4,0.0);	

      h1[ name ]->Fit(f1[name],"R");
      h1[ name ]->SetTitle(yName);
      fit_name[ name ] = yName;
      f1[ name ]->SetTitle(yName);

      f1_vec.push_back( f1[name] );
      name_vec.push_back( name );

    }

}

void ProcessEmbedding::PlotSlices()
{

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  int nSlices = f1.size(); 
  int nPerPage = 9;
  int nPages = f1.size()/nPerPage;
  if ( nPages % nPerPage != 0 ) nPages++;
  
  cout << "nSlices  = " << nSlices << endl;
  cout << "nPerPage = " << nPerPage << endl;
  cout << "nPages   = " << nPages << endl;

  std::vector<TCanvas*> can_vec(nPages,(TCanvas*) NULL);
  int canIndex=0;

  for ( auto &v : can_vec)
    {
      v = new TCanvas(TString::Format("can_%i",canIndex++),"",1200,900);
      v->Divide(3,3,0.0,0.0);
    }

  canIndex=0;
  int padIndex=1;

  for ( int i=2;i<f1_vec.size();i++)
    {

      TF1  * fit  = f1_vec[i];
      TGraphAsymmErrors * hist = h1[ name_vec[i] ];

      can_vec[ canIndex ]->cd( padIndex );
      if ( padIndex == 3 || padIndex == 6 || padIndex == 9 ) gPad->SetRightMargin(0.02);
      //      if ( padIndex > 6 )
      //	{ 
	

      hist->GetYaxis()->SetRangeUser(0.01,1.19);
      hist->GetXaxis()->SetLimits(0.01,2.45);
      

      hist->GetXaxis()->SetNdivisions(508);
      //      hist->GetXaxis()->SetLimits(0.01,2.45);



	  //	}
	  //      else
	  //	{
	  //	  hist->GetYaxis()->SetRangeUser(0.01,1.19);
	  //	  hi/st->GetXaxis()->SetRangeUser(0.01,2.45);
	  //	}

      //      gPad->SetGrid(1,1);
      gPad->SetTicks(1,1);
      hist->SetMarkerStyle( kFullCircle );
      hist->SetMarkerSize( 0.5 );
      hist->GetXaxis()->SetTitle("p_{T} (GeV)");
      hist->GetYaxis()->SetTitle("Eff.");
      //      hist->SetStats(false);
      hist->Draw("ap");
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

  _canvases = can_vec;

}

