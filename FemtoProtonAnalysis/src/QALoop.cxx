#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <numeric>
#include <algorithm>

#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "StFemtoTrack.h"
#include "StFemtoEvent.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3D.h"

#include "InputParameterList.h"
#include "analysisUtil.h"

#include "QALoop.h"


using namespace std;

int SimplePlots(
		TChain * tc,
		long int nentries,
		InputParameterList & pl,
		TString outfilename
		)
{

  map<TString,TH2F*> th2f;
  map<TString,TH3F*> th3f;

  th2f["tpc_eta_pt"]                  = new TH2F("tpc_eta_pt",";#eta;p_{T}",5000,-2.5,0,1000,0,5);
  th2f["tpc_eta_pt_tof"]              = new TH2F("tpc_eta_pt_tof",";#eta;p_{T}",5000,-2.5,0,1000,0,5);
  th2f["tpc_eta_pt_tof_accept"]        = new TH2F("tpc_eta_pt_tof_accept",";#eta;p_{T}",5000,-2.5,0,1000,0,5);
  th2f["tpc_eta_pt_tof_accept_pMom2"] = new TH2F("tpc_eta_pt_tof_accept_pMom2",";#eta;p_{T}",5000,-2.5,0,1000,0,5);

  th3f["fxttof_epd_pipdu"] = new TH3F("fxttof_epd_pipdu","fxttof_epd_pipdu;fxttof;epd;pipdu",200,-0.5,199.5,300,-0.5,299.5,200,-0.5,199.5);
  th3f["fxt3_epd_pipdu"]   = new TH3F("fxt3_epd_pipdu","fxt3_epd_pipdu;fxt3;epd;pipdu",200,-0.5,199.5,300,-0.5,299.5,200,-0.5,199.5);
  th3f["fxt_epd_pipdu"]    = new TH3F("fxt_epd_pipdu","fxt_epd_pipdu;fxt;epd;pipdu",400,-0.5,399.5,300,-0.5,299.5,200,-0.5,199.5);
  th3f["fxt_fxttof_fxt3"]  = new TH3F("fxt_fxttof_fxt3","fxt_fxttof_fxt3;fxt;fxttof;fxt3",400,-0.5,399.5,200,-0.5,199.5,200,-0.5,199.5);
  th3f["fxt_fxttof_fxtpro"]  = new TH3F("fxt_fxttof_fxtpro","fxt_fxttof_fxtpro;fxt;fxttof;fxtpro",400,-0.5,399.5,200,-0.5,199.5,200,-0.5,199.5);

  //  th2f["fxttof_fxt3"]      = new TH2F("fxttof_fxt3","fxttof_fxt3;fxttof;fxt3",200,-0.5,199.5,200,-0.5,199.5);


  th2f["tpc_y_pt"]       = new TH2F("tpc_y_pt",";y;p_{T}",1000,-2.5,2.5,500,0,5);
  th2f["tpc_y_pt_pMom2"] = new TH2F("tpc_y_pt_pMom2",";y;p_{T}",1000,-2.5,2.5,500,0,5);
  //  th2f["tpc_y_pt_pMom2"] = new TH2F("tpc_y_pt_pMom2",";y;p_{T}",1000,-2.5,2.5,500,0,5);
  th2f["tpc_y_pt_tof"]   = new TH2F("tpc_y_pt_tof",";y;p_{T}",1000,-2.5,2.5,500,0,5);

  StFemtoEvent * event = new StFemtoEvent();  
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  //Proton Mass
  double mass=0.938272;

  for (int iEntry=0;iEntry<nentries;iEntry++)
    {
      
      tc->GetEntry(iEntry);

      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
      double fxttof = event->GetFxtMultTofMatch();      

      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      th3f["fxttof_epd_pipdu"]->Fill(fxttof,nmip,pipdu);
      th3f["fxt3_epd_pipdu"]->Fill(fxt3,nmip,pipdu);
      th3f["fxt_epd_pipdu"]->Fill(fxt,nmip,pipdu);

      th3f["fxt_fxttof_fxt3"]->Fill(fxt,fxttof,fxt3);
      th3f["fxt_fxttof_fxtpro"]->Fill(fxt,fxttof,fxt-fxt3);
      //      th2f["fxttof_fxt3"]->Fill(fxttof,fxt3);

      pair<bool,bool> epd_tof = PileUpBadEvent(fxttof,(int)fxt,nmip,pipdu);
      //      if ( !epd_tof.first || !epd_tof.second ) continue;
      
      for (int iTrack=0;iTrack<event->GetEntries();iTrack++)
	{
	  StFemtoTrack trk = event->GetFemtoTrack(iTrack);
	  
	  //Get Track Variables
	  double pt = trk.GetPt();
	  double pz = trk.GetPz();
	  Double_t p      = sqrt( pow(pt,2) + pow(pz,2) );
	  double nsigpro = trk.GetNSigmaProton() - getProtonDedxMeanShift(p);
		
	  double tofmass = trk.GetTofMass();
	  double tofmass2 = (tofmass < 0 ) ?  -1 : tofmass*tofmass;

	  //calculate the energy and rapidity
	  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(mass,2) );

	  Double_t rap = 0.5 * log( ( energy + pz) / (energy - pz) );
	  Double_t eta = 0.5 * log( ( p      + pz) / (p      - pz) );
	  rap = fabs(rap) - 1.04;
	  
	  bool tofmatch=false;
	  bool tofmatchproton=false;
	  if ( tofmass > 0 ) tofmatch=true;
	  if ( tofmass2 > 0.6 && tofmass2 < 1.2) tofmatchproton=true;

	  th2f["tpc_eta_pt"]->Fill( eta, pt);
	  if (tofmatch) th2f["tpc_eta_pt_tof"]->Fill( eta, pt);
	  th2f["tpc_y_pt"]->Fill( rap, pt);
	  if (tofmatchproton || p<2.0) 
	    {
	      if ( fabs(nsigpro) < 3 ) th2f["tpc_y_pt_pMom2"]->Fill( rap, pt);
	    }

	  if (tofmatch) th2f["tpc_y_pt_tof"]->Fill( rap, pt);

	  if ( tofmatch && ( fabs(rap) < 0.5 ) )
	    {
	      th2f["tpc_eta_pt_tof_accept"]->Fill( eta, pt );
	      if ( p > 2.0 ) th2f["tpc_eta_pt_tof_accept_pMom2"]->Fill( eta , pt );
	    }
	}

    }

  TFile * outfile = new TFile(outfilename,"recreate");
  outfile->cd();

  for ( auto &h : th2f ) { h.second->Write(); } 
  for ( auto &h : th3f ) { h.second->Write(); } 


}




int QAplots(
	    TChain * tc,
	    long int nentries,
	    InputParameterList & pl,
	    TString outfilename
	    )
{

  //  map<TString,TH2F*> th2f;
  TFile * outfile = new TFile(outfilename,"recreate");
  //  std::map<TString,TH3F*> th3f;

  /*
  th3f["fxt3_dca_nsigmaproton"] = new TH3F("fxt3_dca_nsigmaproton",";FxtMult3;DCA (cm);n#SigmaProton",400,-0.5,399.5,1000,0,10,10000,-10,10);
  th3f["fxt3_dca_pt"] = new TH3F("fxt3_dca_pt",";FxtMult3;DCA (cm);p_{T}",400,-0.5,399.5,1000,0,10,10000,0,5);
  th3f["fxt3_dca_eta"] = new TH3F("fxt3_dca_eta",";FxtMult3;DCA (cm);#eta",400,-0.5,399.5,1000,0,10,10000,-4,0);
  th3f["fxt3_dca_nhitsfit"] = new TH3F("fxt3_dca_nhitsfit",";FxtMult3;DCA (cm);nHitsFit",400,-0.5,399.5,1000,0,10,100,-0.5,99.5);
  th3f["fxt3_dca_tofmass2"] = new TH3F("fxt3_dca_tofmass2",";FxtMult3;DCA (cm);tofmass2",400,-0.5,399.5,1000,0,10,10000,0,5);
  th3f["fxt3_dca_tofmass2_pcut"] = new TH3F("fxt3_dca_tofmass2_pcut",";FxtMult3;DCA (cm);tofmass2",400,-0.5,399.5,1000,0,10,10000,0,5);
  th3f["fxt3_avgdca_fxt"] = new TH3F("fxt3_avgdca_fxt",";FxtMult3;DCA (cm);FxtMult",400,-0.5,399.5,1000,0,10,400,-0.5,399.5);

  th3f["fxt3_dca_nsigmaproton_aw"] = new TH3F("fxt3_dca_nsigmaproton_aw",";FxtMult3;DCA (cm);n#SigmaProton",400,-0.5,399.5,1000,0,10,10000,-10,10);
  th3f["fxt3_dca_pt_aw"] = new TH3F("fxt3_dca_pt_aw",";FxtMult3;DCA (cm);p_{T}",400,-0.5,399.5,1000,0,10,10000,0,5);
  th3f["fxt3_dca_eta_aw"] = new TH3F("fxt3_dca_eta_aw",";FxtMult3;DCA (cm);#eta",400,-0.5,399.5,1000,0,10,10000,-4,0);
  th3f["fxt3_dca_nhitsfit_aw"] = new TH3F("fxt3_dca_nhitsfit_aw",";FxtMult3;DCA (cm);nHitsFit",400,-0.5,399.5,1000,0,10,100,-0.5,99.5);
  th3f["fxt3_dca_tofmass2_aw"] = new TH3F("fxt3_dca_tofmass2_aw",";FxtMult3;DCA (cm);tofmass2",400,-0.5,399.5,1000,0,10,10000,0,5);
  th3f["fxt3_dca_tofmass2_pcut_aw"] = new TH3F("fxt3_dca_tofmass2_pcut_aw",";FxtMult3;DCA (cm);tofmass2",400,-0.5,399.5,1000,0,10,10000,0,5);
  */

  
  TH1D * fxtmult3_dca3 = new TH1D("fxtmult3","",400,-0.5,399.5);
  TH1D * fxtmult3_dca2p5 = new TH1D("fxtmult3_dca2p5","",400,-0.5,399.5);
  TH1D * fxtmult3_dca2 = new TH1D("fxtmult3_dca2","",400,-0.5,399.5);
  TH1D * fxtmult3_dca1 = new TH1D("fxtmult3_dca1","",400,-0.5,399.5);

  TH1D * fxtmult_dca3 = new TH1D("fxtmult_dca3","",400,-0.5,399.5);
  TH1D * fxtmult_dca2p5 = new TH1D("fxtmult_dca2p5","",400,-0.5,399.5);
  TH1D * fxtmult_dca2 = new TH1D("fxtmult_dca2","",400,-0.5,399.5);
  TH1D * fxtmult_dca1 = new TH1D("fxtmult_dca1","",400,-0.5,399.5);

  TH3D * th3f_cent_dca_nsigmaproton = new TH3D("cent_dca_nsigmaproton",";FxtMult3;DCA (cm);n#SigmaProton",9,-0.5,8.5,300,0,3,1000,-10,10);
  TH3D * th3f_cent_dca_pt = new TH3D("cent_dca_pt",";FxtMult3;DCA (cm);p_{T}",9,-0.5,8.5,300,0,3,1000,0,5);
  TH3D * th3f_cent_dca_eta = new TH3D("cent_dca_eta",";FxtMult3;DCA (cm);#eta",9,-0.5,8.5,300,0,3,1000,-4,0);
  TH3D * th3f_cent_dca_nhitsfit = new TH3D("cent_dca_nhitsfit",";FxtMult3;DCA (cm);nHitsFit",9,-0.5,8.5,300,0,3,100,-0.5,99.5);
  TH3D * th3f_cent_dca_tofmass2 = new TH3D("cent_dca_tofmass2",";FxtMult3;DCA (cm);tofmass2",9,-0.5,8.5,300,0,3,1000,0,5);
  TH3D * th3f_cent_dca_tofmass2_pcut = new TH3D("cent_dca_tofmass2_pcut",";FxtMult3;DCA (cm);tofmass2",9,-0.5,8.5,300,0,3,1000,0,5);
  TH3D * th3f_fxt3_avgdca_fxt = new TH3D("fxt3_avgdca_fxt",";FxtMult3;DCA (cm);FxtMult",400,-0.5,399.5,300,0,3,400,-0.5,399.5);

  TH3D * th3f_cent_dca_nsigmaproton_aw = new TH3D("cent_dca_nsigmaproton_aw",";FxtMult3;DCA (cm);n#SigmaProton",9,-0.5,8.5,300,0,3,1000,-10,10);
  TH3D * th3f_cent_dca_pt_aw = new TH3D("cent_dca_pt_aw",";FxtMult3;DCA (cm);p_{T}",9,-0.5,8.5,300,0,3,1000,0,5);
  TH3D * th3f_cent_dca_eta_aw = new TH3D("cent_dca_eta_aw",";FxtMult3;DCA (cm);#eta",9,-0.5,8.5,300,0,3,1000,-4,0);
  TH3D * th3f_cent_dca_nhitsfit_aw = new TH3D("cent_dca_nhitsfit_aw",";FxtMult3;DCA (cm);nHitsFit",9,-0.5,8.5,300,0,3,100,-0.5,99.5);
  TH3D * th3f_cent_dca_tofmass2_aw = new TH3D("cent_dca_tofmass2_aw",";FxtMult3;DCA (cm);tofmass2",9,-0.5,8.5,300,0,3,1000,0,5);
  TH3D * th3f_cent_dca_tofmass2_pcut_aw = new TH3D("cent_dca_tofmass2_pcut_aw",";FxtMult3;DCA (cm);tofmass2",9,-0.5,8.5,300,0,3,1000,0,5);
  
  StFemtoEvent * event = new StFemtoEvent();  
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  //Proton Mass
  double mass=0.938272;

  for (int iEntry=0;iEntry<nentries;iEntry++)
    {
      
      tc->GetEntry(iEntry);
      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
      double fxttof = event->GetFxtMultTofMatch();      


      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      double fxt_dca3=0;
      double fxt_dca2p5=0;
      double fxt_dca2=0;
      double fxt_dca1=0;

      double fxt3_dca3=0;
      double fxt3_dca2p5=0;
      double fxt3_dca2=0;
      double fxt3_dca1=0;

      int centBin = CentBin3(fxt3);
      if ( centBin == -1 ) centBin=8;
      std::vector<double> avg_dca;

      for (int iTrack=0;iTrack<event->GetEntries();iTrack++)
	{
	  StFemtoTrack trk = event->GetFemtoTrack(iTrack);
	  

	  //Get Track Variables
	  short nhitsfit = trk.GetNHitsFit();
	  double dcaZ  = trk.GetDcaZ();
	  double dcaX  = trk.GetDcaX();
	  double dcaY  = trk.GetDcaY();
	  double dca = sqrt(pow(dcaX,2) + pow(dcaY,2) + pow(dcaZ,2));

	  if ( dca < 3 ) fxt_dca3+=1;
	  if ( dca < 2.5 ) fxt_dca2p5+=1;
	  if ( dca < 2 ) fxt_dca2+=1;
	  if ( dca < 1 ) fxt_dca1+=1;

	  double pt = trk.GetPt();
	  double pz = trk.GetPz();
	  Double_t p      = sqrt( pow(pt,2) + pow(pz,2) );
	  double nsigpro = trk.GetNSigmaProton() - getProtonDedxMeanShift(p);

	  if ( nsigpro < -3 )
	    {
	      if ( dca < 3 ) fxt3_dca3+=1;
	      if ( dca < 2.5 ) fxt3_dca2p5+=1;
	      if ( dca < 2 ) fxt3_dca2+=1;
	      if ( dca < 1 ) fxt3_dca1+=1;
	    }
		
	  double tofmass = trk.GetTofMass();
	  double tofmass2 = (tofmass < 0 ) ?  -1 : tofmass*tofmass;

	  //calculate the energy and rapidity
	  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(mass,2) );

	  Double_t rap = 0.5 * log( ( energy + pz) / (energy - pz) );
	  Double_t eta = 0.5 * log( ( p      + pz) / (p      - pz) );
	  rap = fabs(rap) - 1.049;
	  
	  bool tofmatch=false;
	  bool tofmatchproton=false;
	  bool goodProton=false;
	  if ( tofmass > 0 ) tofmatch=true;
	  if ( tofmass2 > 0.6 && tofmass2 < 1.2) tofmatchproton=true;


	  if (tofmatchproton || p<2.0) 
	    {
	      if ( fabs(nsigpro) < 3 )
		{
		  if ( rap > -0.5 && rap < 0 )
		    {
		      if (pt > 0.4 && pt < 2.0 )
			{
			  goodProton=true;
			}
		    }	  
		}
	    }

	  th3f_cent_dca_nsigmaproton->Fill( centBin, dca, nsigpro);
	  th3f_cent_dca_pt->Fill(centBin, dca, pt);
	  th3f_cent_dca_eta->Fill( centBin, dca, eta);
	  th3f_cent_dca_nhitsfit->Fill( centBin, dca, nhitsfit);
	  th3f_cent_dca_tofmass2->Fill( centBin, dca, tofmass2);
	  if ( p > 2.0 ) th3f_cent_dca_tofmass2_pcut->Fill( centBin, dca, tofmass2);

	  if ( goodProton )
	    {
	      th3f_cent_dca_nsigmaproton_aw->Fill( centBin, dca, nsigpro);
	      th3f_cent_dca_pt_aw->Fill( centBin, dca, pt);
	      th3f_cent_dca_eta_aw->Fill( centBin, dca, eta);
	      th3f_cent_dca_nhitsfit_aw->Fill( centBin, dca, nhitsfit);
	      th3f_cent_dca_tofmass2_aw->Fill( centBin, dca, tofmass2);
	      if ( p > 2.0 ) th3f_cent_dca_tofmass2_pcut_aw->Fill( centBin, dca, tofmass2);
	    }
	  avg_dca.push_back(dca);
	
	}

      double avgDca = accumulate(avg_dca.begin(),avg_dca.end(),0.0);
      th3f_fxt3_avgdca_fxt->Fill( fxt3,avgDca,fxt_dca3);

      fxtmult_dca3->Fill(fxt_dca3);
      fxtmult_dca2p5->Fill(fxt_dca2p5);
      fxtmult_dca2->Fill(fxt_dca2);
      fxtmult_dca1->Fill(fxt_dca1);

      fxtmult3_dca3->Fill(fxt3_dca3);
      fxtmult3_dca2p5->Fill(fxt3_dca2p5);
      fxtmult3_dca2->Fill(fxt3_dca2);
      fxtmult3_dca1->Fill(fxt3_dca1);
  
    }

  outfile->cd();

  fxtmult_dca3->Write();
  fxtmult_dca2p5->Write();
  fxtmult_dca2->Write();
  fxtmult_dca1->Write();

  fxtmult3_dca3->Write();
  fxtmult3_dca2p5->Write();
  fxtmult3_dca2->Write();
  fxtmult3_dca1->Write();

  th3f_cent_dca_nsigmaproton->Write();
  th3f_cent_dca_pt->Write();
  th3f_cent_dca_eta->Write();
  th3f_cent_dca_nhitsfit->Write();
  th3f_cent_dca_tofmass2->Write();
  th3f_cent_dca_tofmass2_pcut->Write();
  th3f_fxt3_avgdca_fxt->Write();

  th3f_cent_dca_nsigmaproton_aw->Write();
  th3f_cent_dca_pt_aw->Write();
  th3f_cent_dca_eta_aw->Write();
  th3f_cent_dca_nhitsfit_aw->Write();
  th3f_cent_dca_tofmass2_aw->Write();
  th3f_cent_dca_tofmass2_pcut_aw->Write();

  //  for ( auto &h : th2f ) { h.second->Write(); } 
  //  for ( auto &h : th3f ) { h.second->Write(); } 


}


int MakeMultHists(
		  TChain * tc,
		  long int nentries,
		  InputParameterList & pl,
		  TString outfilename
		  )
{


  TFile * outfile = new TFile(outfilename,"recreate");

  TH1D * fxtmult3 = new TH1D("fxtmult3","",400,-0.5,399.5);
  TH1D * fxtmult = new TH1D("fxtmult","",400,-0.5,399.5);

  StFemtoEvent * event = new StFemtoEvent();  
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  //Proton Mass

  for (int iEntry=0;iEntry<nentries;iEntry++)
    {
      
      tc->GetEntry(iEntry);
      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
      double fxttof = event->GetFxtMultTofMatch();      

      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      fxtmult->Fill(fxt);
      fxtmult3->Fill(fxt3);
      
    }

  outfile->cd();

  fxtmult->Write();
  fxtmult3->Write();

}


int QALoop(
	   TChain * tc,
	   long int nentries,
	   InputParameterList & pl,
	   TFile * outfile
	   )
{

  map<TString,TH1F*> th1f;
  map<TString,TH2F*> th2f;
  map<TString,TProfile*> prof;
  map<TString,TProfile2D*> prof2d;
  std::vector<int> runNumVec = {19153029, 19153031, 19153033, 19153034, 19153035, 19153036, 19153037, 19153042, 19153043,
				19153044, 19153050, 19153051, 19153052, 19153053, 19153054, 19153055, 19153056, 19153057,
				19153058, 19153059, 19153061, 19153062, 19153063, 19153064, 19153066, 19154001, 19154002,
				19154005, 19154007, 19154027, 19154028, 19154029, 19154030, 19154031, 19154032, 19154036,
				19154037, 19154038, 19154039, 19154040, 19154041, 19154044, 19154045, 19154046, 19154047,
				19154048, 19154049, 19154052, 19154053, 19154054, 19154055, 19154056, 19154057, 19154058,
				19154061, 19154063, 19154064, 19154065, 19154066, 19154067, 19155001, 19155003, 19155004,
				19155005, 19155006, 19155008, 19155009, 19155010, 19155011, 19155016, 19155017, 19155018,
				19155019, 19155020, 19155021, 19155022};
  
  //Event TH1F
  th1f["Vxy"] = new TH1F("Vxy",";Vr from (-2,0)cm",500,-5,5);
  th1f["Vz"] = new TH1F("Vz",";Vz (cm)",1000,195,205);

  th1f["refmult"] = new TH1F("refmult",";FxtMult",300,0,300);
  th1f["refmult3"] = new TH1F("refmult3",";FxtMult3",300,0,300);
  th1f["refmult_pro"] = new TH1F("refmult_pro",";FxtMult Protons",300,0,300);

  //Event TH1F no cuts
  th1f["Vxy_nc"] = new TH1F("Vxy_nc",";Vr from (-2,0)cm",500,-5,5);
  th1f["Vz_nc"] = new TH1F("Vz_nc",";Vz (cm)",1000,195,205);

  th1f["refmult_nc"] = new TH1F("refmult_nc",";FxtMult",300,0,300);
  th1f["refmult3_nc"] = new TH1F("refmult3_nc",";FxtMult3",300,0,300);

  //Event TH2F
  th2f["Vxy_Vz"] = new TH2F("VxyVz",";Vr (cm);Vz (cm)",500,-5,5,500,195,205);
  th2f["refmult_nmip"] = new TH2F("refmult_mnmip",";FxtMult;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["refmult3_nmip"] = new TH2F("refmult3_mnmip",";FxtMult3;EPD #Sigma nMip",300,0,300,200,0,200);

  th2f["refmult_pipdu"] = new TH2F("refmult_pipdu",";FxtMult;TOF #pi + p + du",300,0,300,300,0,300);
  th2f["refmult3_pipdu"] = new TH2F("refmult3_pipdu",";FxtMult3;TOF #pi + p + du",300,0,300,300,0,300);

  th2f["refmult_refmult3"] = new TH2F("refmult_refmult3",";FxtMult;FxtMult3",300,0,300,300,0,300);

  //Event TH2F no cuts
  th2f["Vxy_Vz_nc"] = new TH2F("VxyVz_nc",";Vr (cm);Vz (cm)",500,-5,5,500,195,205);
  th2f["refmult_nmip_nc"] = new TH2F("refmult_mnmip_nc",";FxtMult;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["refmult3_nmip_nc"] = new TH2F("refmult3_mnmip_nc",";FxtMult3;EPD #Sigma nMip",300,0,300,200,0,200);

  th2f["refmult_pipdu_nc"] = new TH2F("refmult_pipdu_nc",";FxtMult;TOF #pi + p + du",300,0,300,300,0,300);
  th2f["refmult3_pipdu_nc"] = new TH2F("refmult3_pipdu_nc",";FxtMult3;TOF #pi + p + du",300,0,300,300,0,300);

  th2f["refmult_refmult3_nc"] = new TH2F("refmult_refmult3_nc",";FxtMult;FxtMult3",300,0,300,300,0,300);
  
  //Track TH1F
  th1f["pt"] = new TH1F("pt",";p_{T} (GeV/c)",500,0,5);
  th1f["pz"] = new TH1F("pz",";p_{Z} (GeV/c)",500,-10,10);
  th1f["nhitsfit"] = new TH1F("nhitsfit",";n HitsFit",80,0,80);
  th1f["dcaz"] = new TH1F("dcaz",";<dcaZ> (cm)",100,-5,5);
  th1f["dcaxy"] = new TH1F("dcaxy",";<dcaXY> (cm)",100,-5,5);
  th1f["dcaxypro"] = new TH1F("dcaxypro",";<dcaXY> (cm)",100,-5,5);
  th1f["dcaxy2"] = new TH1F("dcaxy2",";<dcaXY> (cm)",100,-5,5);
  th1f["nsigproton"] = new TH1F("nsigproton",";N #sigma Proton",100,-5,5);
  th1f["tofmass2"] = new TH1F("tofmass2",";m^{2} (GeV/c^{2})^{2}",500,0,5);

  //Track TH1F no cuts
  th1f["pt_nc"] = new TH1F("pt_nc",";p_{T} (GeV/c)",500,0,5);
  th1f["pz_nc"] = new TH1F("pz_nc",";p_{T} (GeV/c)",500,-10,10);
  th1f["nhitsfit_nc"] = new TH1F("nhitsfit_nc",";n HitsFit",80,0,80);
  th1f["dcaz_nc"] = new TH1F("dcaz_nc",";dcaZ (cm)",100,-5,5);
  th1f["dcaxy_nc"] = new TH1F("dcaxy_nc",";<dcaXY> (cm)",100,-5,5);
  th1f["nsigproton_nc"] = new TH1F("nsigproton_nc",";N #sigma Proton",100,-5,5);
  th1f["tofmass2_nc"] = new TH1F("tofmass2_nc",";m^{2} (GeV/c^{2})^{2}",500,0,5);

  //Track TH2F
  th2f["tofmass2_nsigproton"] = new TH2F("tofmass2_nsigproton",";m^{2} (GeV/c^{2})^2; N #sigma Proton",500,0,5,100,-5,5);
  th2f["y_pt_pos"] = new TH2F("y_pt_pos",";y-y_{cm}; p_{T} (GeV/c)",500,-1.5,1.5,500,0,5);
  th2f["y_pt_neg"] = new TH2F("y_pt_neg",";y-y_{cm}; p_{T} (GeV/c)",500,-1.5,1.5,500,0,5);

  th2f["dcaxy_npro"] = new TH2F("dcaxy_npro",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxy_npro_epdcut"] = new TH2F("dcaxy_npro_epdcut",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxy_npro_tofcut"] = new TH2F("dcaxy_npro_tofcut",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxy_npro_all"] = new TH2F("dcaxy_npro_all",";<dcaXY> (cm)",1000,-5,5,200,0,200);

  th2f["dcaxypro_npro"] = new TH2F("dcaxypro_npro",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_epdcut"] = new TH2F("dcaxypro_npro_epdcut",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_tofcut"] = new TH2F("dcaxypro_npro_tofcut",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_all"] = new TH2F("dcaxypro_npro_all",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  
  th2f["dcaxy_fxt3"] = new TH2F("dcaxy_fxt3",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_epdcut"] = new TH2F("dcaxy_fxt3_epdcut",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_tofcut"] = new TH2F("dcaxy_fxt3_tofcut",";<dcaXY> (cm)",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_all"] = new TH2F("dcaxy_fxt3_all",";<dcaXY> (cm)",1000,-5,5,200,0,200);
 
  th2f["dcaxypro_npro_cent5"] = new TH2F("dcaxypro_npro_cent5",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_all_cent5"] = new TH2F("dcaxypro_npro_all_cent5",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200); 
  th2f["dcaxy_fxt3_cent5"] = new TH2F("dcaxy_fxt3_cent5",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_all_cent5"] = new TH2F("dcaxy_fxt3_all_cent5",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);

  th2f["dcaxypro_npro_cent10"] = new TH2F("dcaxypro_npro_cent10",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_all_cent10"] = new TH2F("dcaxypro_npro_all_cent10",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_cent10"] = new TH2F("dcaxy_fxt3_cent10",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_all_cent10"] = new TH2F("dcaxy_fxt3_all_cent10",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);

  th2f["dcaxypro_npro_cent20"] = new TH2F("dcaxypro_npro_cent20",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_all_cent20"] = new TH2F("dcaxypro_npro_all_cent20",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);  
  th2f["dcaxy_fxt3_cent20"] = new TH2F("dcaxy_fxt3_cent20",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_all_cent20"] = new TH2F("dcaxy_fxt3_all_cent20",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);

  th2f["dcaxypro_npro_cent50"] = new TH2F("dcaxypro_npro_cent50",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);
  th2f["dcaxypro_npro_all_cent50"] = new TH2F("dcaxypro_npro_all_cent50",";<dcaXY> (cm);nProtons",1000,-5,5,200,0,200);  
  th2f["dcaxy_fxt3_cent50"] = new TH2F("dcaxy_fxt3_cent50",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);
  th2f["dcaxy_fxt3_all_cent50"] = new TH2F("dcaxy_fxt3_all_cent50",";<dcaXY> (cm);FxtMult3",1000,-5,5,200,0,200);

  
  //Track TH2F no cuts
  th2f["tofmass2_nsigproton_nc"] = new TH2F("tofmass2_nsigproton_nc",";m^{2} (GeV/c^{2})^2; N #sigma Proton",500,0,5,100,-5,5);
  th2f["y_pt_pos_nc"] = new TH2F("y_pt_pos_nc",";y-y_{cm}; p_{T} (GeV/c)",500,-1.5,1.5,500,0,5);
  th2f["tpc"] = new TH2F("tpc",";y-y_{cm}; p_{T} (GeV/c)",500,-1.5,1.5,500,0,5);
  th2f["tpc_tofmatch"] = new TH2F("tpc_tofmatch",";y-y_{cm}; p_{T} (GeV/c)",500,-1.5,1.5,500,0,5);
  
  th2f["y_pt_neg_nc"] = new TH2F("y_pt_neg_nc",";y-y_{cm}; p_{T} (GeV/c)",500,-1.5,1.5,500,0,5);

  //Profiles
  prof["Vxy_prof"] = new TProfile("Vxy_prof",";Run Index;Vr from (-2,0)cm",500,-0.5,499.5);
  prof["Vz_prof"] = new TProfile("Vz_prof",";Run Index;Vz (cm)",500,-0.5,499.5);

  prof["refmult_prof"] = new TProfile("refmult_prof",";Run Index;FxtMult",500,-0.5,499.5);
  prof["refmult3_prof"] = new TProfile("refmult3_prof",";Run Index;FxtMult3",500,-0.5,499.5);

  prof["pt_prof"] = new TProfile("pt_prof",";Run Index;<p_{T}> (GeV/c)",500,-0.5,499.5);
  prof["pz_prof"] = new TProfile("pz_prof",";Run Index;<p_{Z}> (GeV/c)",500,-0.5,499.5);
  prof["nhitsfit_prof"] = new TProfile("nhitsfit_prof",";Run Index;<n HitsFit>",500,-0.5,499.5);
  prof["dcaz_prof"] = new TProfile("dcaz_prof",";Run Index;<dcaZ> (cm)",500,-0.5,499.5);
  prof["dcaxy_prof"] = new TProfile("dcaxy_prof",";Run Index;<dcaXY> (cm)",500,-0.5,499.5);

  //Profiles 2D
  prof2d["dcaxy_fxt3_prof"] = new TProfile2D("dcaxy_fxt3_prof",";Run Index;FxtMult3",500,-0.5,499.5,200,0,200);
  prof2d["dcaxy_fxt3_all_prof"] = new TProfile2D("dcaxy_fxt3_all_prof",";Run Index;FxtMult3",500,-0.5,499.5,200,0,200);
  prof2d["dcaxy_fxt3_epdcut_prof"] = new TProfile2D("dcaxy_fxt3_epdcut_prof",";Run Index;FxtMult3",500,-0.5,499.5,200,0,200);
  prof2d["dcaxy_fxt3_tofcut_prof"] = new TProfile2D("dcaxy_fxt3_tofcut_prof",";Run Index;FxtMult3",500,-0.5,499.5,200,0,200);

  prof2d["dcaxy_npro_prof"] = new TProfile2D("dcaxy_npro_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  prof2d["dcaxy_npro_all_prof"] = new TProfile2D("dcaxy_npro_all_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  prof2d["dcaxy_npro_epdcut_prof"] = new TProfile2D("dcaxy_npro_epdcut_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  prof2d["dcaxy_npro_tofcut_prof"] = new TProfile2D("dcaxy_npro_tofcut_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);

  prof2d["dcaxypro_npro_prof"] = new TProfile2D("dcaxypro_npro_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  prof2d["dcaxypro_npro_all_prof"] = new TProfile2D("dcaxypro_npro_all_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  prof2d["dcaxypro_npro_epdcut_prof"] = new TProfile2D("dcaxypro_npro_epdcut_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  prof2d["dcaxypro_npro_tofcut_prof"] = new TProfile2D("dcaxypro_npro_tofcut_prof",";Run Index;nProtons",500,-0.5,499.5,200,0,200);
  
  StFemtoEvent * event = new StFemtoEvent();  
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  //Proton Mass
  double mass=0.938272;

  long int onePercent = nentries/100;
  int percent = 1;

  clock_t t0 = clock();

  
  cout << "Total Number of entries =" << nentries << endl; 
  
  for (int iEntry=0;iEntry<nentries;iEntry++)
    {

      if ( iEntry % onePercent == 0 && iEntry != 0 )
	{
	  clock_t t = clock() - t0;
	  cout << percent  << "%" << " iEntry =" << iEntry << endl;
	  cout << "                          Run Time = " << ((float)t/CLOCKS_PER_SEC) << " secs" << endl;
	  cout << "          Estimated Time Remaining = " << ((float)t/CLOCKS_PER_SEC)*(100 - percent)/percent << endl;
	  percent = percent + 1;
	}
      
      tc->GetEntry(iEntry);

      int runIndex = GetRunIndex( runNumVec , event->GetRunID() );
      
      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
      //Fill no cut Event histograms

      th1f["Vxy_nc"]->Fill(vr);
      th1f["Vz_nc"]->Fill(vz);

      th2f["refmult_refmult3_nc"]->Fill(fxt,fxt3);
      th1f["refmult_nc"]->Fill(fxt);
      th1f["refmult3_nc"]->Fill(fxt3);
      
      th2f["Vxy_Vz_nc"]->Fill(vr,vz);
      th2f["refmult_nmip_nc"]->Fill(fxt,nmip);
      th2f["refmult3_nmip_nc"]->Fill(fxt3,nmip);
      
      th2f["refmult_pipdu_nc"]->Fill(fxt,pipdu);
      th2f["refmult3_pipdu_nc"]->Fill(fxt3,pipdu);

      prof["Vxy_prof"]->Fill(runIndex, vr);
      prof["Vz_prof"]->Fill(runIndex, vz);
      
      prof["refmult_prof"]->Fill(runIndex, fxt);
      prof["refmult3_prof"]->Fill(runIndex, fxt3);
      
      //Make Event Cuts
      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      
      pair<bool,bool> epd_tof = PileUpBadEvent((int)fxt,nmip,pipdu);

      /*
      if (epd_tof.first)
	{
	  th2f["refmult_nmip"]->Fill(fxt,nmip);
	  th2f["refmult3_nmip"]->Fill(fxt3,nmip);
	}
      if (epd_tof.second)
	{
	  th2f["refmult_pipdu"]->Fill(fxt,pipdu);
	  th2f["refmult3_pipdu"]->Fill(fxt3,pipdu);
	}
			      */
      //      if (!epd_tof.first || !epd_tof.second) continue;

      th1f["Vxy"]->Fill(vr);
      th1f["Vz"]->Fill(vz);
      th2f["Vxy_Vz"]->Fill(vr,vz);
      th2f["refmult_refmult3"]->Fill(fxt,fxt3);
      th1f["refmult"]->Fill(fxt);
      th1f["refmult3"]->Fill(fxt3);
      th1f["refmult_pro"]->Fill(fxt-fxt3);

      std::vector<double> avg_dcaXY;
      std::vector<double> avg_dcaXYpro;
      std::vector<double> avg_dcaXY2;
      std::vector<double> avg_dcaZ;
      int ntot=0;
      int npro=0;
      
      for (int iTrack=0;iTrack<event->GetEntries();iTrack++)
	{
	  StFemtoTrack trk = event->GetFemtoTrack(iTrack);
	  
	  //Get Track Variables
	  double pt = trk.GetPt();
	  double pz = trk.GetPz();
	  short nhitsfit = trk.GetNHitsFit();
	  short charge = trk.GetCharge();
	  double dcaZ  = trk.GetDcaZ();
	  double dcaX  = trk.GetDcaX();
	  double dcaY  = trk.GetDcaY();
	  double nsigpro = trk.GetNSigmaProton();
	  double nsigpi = trk.GetNSigmaPion();
	  double tofmass = trk.GetTofMass();

	  //calculate the energy and rapidity
	  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(mass,2) );
	  Double_t rap = 0.5 * log( ( energy + pz) / (energy - pz) );
	  rap = fabs(rap);

	  prof["pt_prof"]->Fill( runIndex, pt);
	  prof["pz_prof"]->Fill( runIndex, pz);
	  prof["nhitsfit_prof"]->Fill( runIndex, nhitsfit);
	  
	  th1f["pt_nc"]->Fill(pt);
	  th1f["pz_nc"]->Fill(pz);
	  th1f["nhitsfit_nc"]->Fill(nhitsfit);
	  //	  th1f["dcaz_nc"]->Fill(dcaZ);
	  //	  th1f["dcaxy_nc"]->Fill(dcaXY);
	  th1f["nsigproton_nc"]->Fill(nsigpro);
	  th1f["tofmass2_nc"]->Fill(tofmass*tofmass);

	  th2f["tofmass2_nsigproton_nc"]->Fill(tofmass*tofmass,nsigpro);
	  if (charge >0 ) th2f["y_pt_pos_nc"]->Fill(rap-1.04,pt);
	  if (charge <0 ) th2f["y_pt_neg_nc"]->Fill(rap-1.04,pt);

	  if ( fabs(nsigpro) < 3 ) npro += 1;
	  

	  TVector3 pMom = TVector3(trk.GetPx(),trk.GetPy(),trk.GetPz()).Unit();
	  float cos1 = pMom.Perp();
	  float dcaD = dcaX*pMom.Y()/cos1  - dcaY*pMom.X()/cos1;
	  float dcaD1 = dcaX*pMom.X()/cos1 + dcaY*pMom.Y()/cos1; 

	  int sign = (dcaD > 0 ) ? 1 : -1;
	  float dcaXY = sign*sqrt(dcaX*dcaX + dcaY*dcaY);
	  //	  cout << "diff = " << dcaD - dcaD1 << endl;
	  
	  avg_dcaXY.push_back( dcaXY );
	  if ( fabs(nsigpro) < 2 ) avg_dcaXYpro.push_back( dcaXY );
	  avg_dcaXY2.push_back( dcaD1 );
	  avg_dcaZ.push_back(dcaZ);
	  ntot++;
	  
	  //Make Track Cuts
	  if ( fabs(nsigpro) > pl.NSigmaProtonCut()) continue;
	  if ( fabs(nsigpi) < pl.NSigmaPionCut()) continue;
	  if ( nhitsfit <= pl.NHitsFitMin()) continue; 

	  if (charge >0 ) th2f["tpc"]->Fill(rap-1.04,pt);
	  if (charge >0 && tofmass*tofmass > 0.6 &&  tofmass*tofmass < 1.2) th2f["tpc_tofmatch"]->Fill(rap-1.04,pt);
	  
	  if ( pt < pl.PtLow() ) continue;
	  if ( pt > pl.PtHigh() ) continue;
	  if ( rap < pl.RapidityLow() ) continue;
	  if ( rap > pl.RapidityHigh() ) continue;

	  th1f["pt"]->Fill(pt);
	  th1f["pz"]->Fill(pz);
	  th1f["nhitsfit"]->Fill(nhitsfit);
	  //	  th1f["dcaz"]->Fill(dcaZ);
	  //	  th1f["dcaxy"]->Fill(dcaXY);
	  th1f["nsigproton"]->Fill(nsigpro);
	  th1f["tofmass2"]->Fill(tofmass*tofmass);

	  th2f["tofmass2_nsigproton"]->Fill(tofmass*tofmass,nsigpro);
	  if (charge >0 ) th2f["y_pt_pos"]->Fill(rap-1.04,pt);
	  if (charge <0 ) th2f["y_pt_neg"]->Fill(rap-1.04,pt);
	  	  
	}

      double avgDcaXY = accumulate(avg_dcaXY.begin(),avg_dcaXY.end(),0.0);
      double avgDcaXYpro = accumulate(avg_dcaXYpro.begin(),avg_dcaXYpro.end(),0.0);
      double avgDcaXY2 = accumulate(avg_dcaXY2.begin(),avg_dcaXY2.end(),0.0);
      double avgDcaZ = accumulate(avg_dcaZ.begin(),avg_dcaZ.end(),0.0);

      if ( ntot == 0 ) continue;
      
      avgDcaXY=avgDcaXY/ntot;
      avgDcaXYpro=avgDcaXY/((float)avg_dcaXYpro.size());
      avgDcaZ=avgDcaZ/ntot;
      
      
      th1f["dcaz"]->Fill(avgDcaZ);
      th1f["dcaxy"]->Fill(avgDcaXY);
      th1f["dcaxypro"]->Fill(avgDcaXYpro);
      th1f["dcaxy2"]->Fill(avgDcaXY2);

      //      cout << "Run Index = " << runIndex << " avgDcaXY = " << avgDcaXY << endl;
      
      if ( isfinite(avgDcaZ) ) prof["dcaz_prof"]->Fill( runIndex, avgDcaZ);
      if ( isfinite(avgDcaXY) )prof["dcaxy_prof"]->Fill( runIndex, avgDcaXY);

      
      th2f["dcaxy_npro"]->Fill(avgDcaXY,npro);
      if (epd_tof.first)                   th2f["dcaxy_npro_epdcut"]->Fill(avgDcaXY,npro);
      if (epd_tof.second)                  th2f["dcaxy_npro_tofcut"]->Fill(avgDcaXY,npro);
      if (epd_tof.first && epd_tof.second) th2f["dcaxy_npro_all"]->Fill(avgDcaXY,npro);

      th2f["dcaxypro_npro"]->Fill(avgDcaXYpro,npro);
      if (epd_tof.first)                   th2f["dcaxypro_npro_epdcut"]->Fill(avgDcaXYpro,npro);
      if (epd_tof.second)                  th2f["dcaxypro_npro_tofcut"]->Fill(avgDcaXYpro,npro);
      if (epd_tof.first && epd_tof.second) th2f["dcaxypro_npro_all"]->Fill(avgDcaXYpro,npro);

      th2f["dcaxy_fxt3"]->Fill(avgDcaXY,fxt3);
      if (epd_tof.first)                   th2f["dcaxy_fxt3_epdcut"]->Fill(avgDcaXY,fxt3);
      if (epd_tof.second)                  th2f["dcaxy_fxt3_tofcut"]->Fill(avgDcaXY,fxt3);
      if (epd_tof.first && epd_tof.second) th2f["dcaxy_fxt3_all"]->Fill(avgDcaXY,fxt3);


      int centBin3 = CentBin3(fxt3);

      
      if (centBin3 == 0 )
	{
	  th2f["dcaxypro_npro_cent5"]->Fill(avgDcaXYpro,npro);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxypro_npro_all_cent5"]->Fill(avgDcaXYpro,npro);

	  th2f["dcaxy_fxt3_cent5"]->Fill(avgDcaXY,fxt3);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxy_fxt3_all_cent5"]->Fill(avgDcaXY,fxt3);	  
	}

      if (centBin3 == 1 )
	{
	  th2f["dcaxypro_npro_cent10"]->Fill(avgDcaXYpro,npro);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxypro_npro_all_cent10"]->Fill(avgDcaXYpro,npro);

	  th2f["dcaxy_fxt3_cent10"]->Fill(avgDcaXY,fxt3);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxy_fxt3_all_cent10"]->Fill(avgDcaXY,fxt3);	  
	}
      
      if (centBin3 == 2 || centBin3 == 3 )
	{
	  th2f["dcaxypro_npro_cent20"]->Fill(avgDcaXYpro,npro);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxypro_npro_all_cent20"]->Fill(avgDcaXYpro,npro);
	  
	  th2f["dcaxy_fxt3_cent20"]->Fill(avgDcaXY,fxt3);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxy_fxt3_all_cent20"]->Fill(avgDcaXY,fxt3);	  
	}
      
      if (centBin3 > 3 && centBin3 <= 9 )
	{
	  th2f["dcaxypro_npro_cent50"]->Fill(avgDcaXYpro,npro);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxypro_npro_all_cent50"]->Fill(avgDcaXYpro,npro);

	  th2f["dcaxy_fxt3_cent50"]->Fill(avgDcaXY,fxt3);
	  if (epd_tof.first && epd_tof.second) th2f["dcaxy_fxt3_all_cent50"]->Fill(avgDcaXY,fxt3);	  
	}


      
      prof2d["dcaxy_fxt3_prof"]->Fill( runIndex, fxt3, avgDcaXY, 1);
      if ( epd_tof.first && epd_tof.second ) prof2d["dcaxy_fxt3_all_prof"]->Fill( runIndex, fxt3, avgDcaXY, 1);
      if ( epd_tof.first )                   prof2d["dcaxy_fxt3_epdcut_prof"]->Fill( runIndex, fxt3, avgDcaXY, 1);
      if ( epd_tof.second )                  prof2d["dcaxy_fxt3_tofcut_prof"]->Fill( runIndex, fxt3, avgDcaXY, 1);
      
      prof2d["dcaxy_npro_prof"]->Fill( runIndex, npro, avgDcaXY, 1);
      if ( epd_tof.first && epd_tof.second ) prof2d["dcaxy_npro_all_prof"]->Fill( runIndex, npro, avgDcaXY, 1);
      if ( epd_tof.first )                   prof2d["dcaxy_npro_epdcut_prof"]->Fill( runIndex, npro, avgDcaXY, 1);
      if ( epd_tof.second )                  prof2d["dcaxy_npro_tofcut_prof"]->Fill( runIndex, npro, avgDcaXY, 1);
      
      prof2d["dcaxypro_npro_prof"]->Fill( runIndex, npro, avgDcaXYpro, 1);
      if ( epd_tof.first && epd_tof.second ) prof2d["dcaxypro_npro_all_prof"]->Fill( runIndex, npro, avgDcaXYpro, 1);
      if ( epd_tof.first )                   prof2d["dcaxypro_npro_epdcut_prof"]->Fill( runIndex, npro, avgDcaXYpro, 1);
      if ( epd_tof.second )                  prof2d["dcaxypro_npro_tofcut_prof"]->Fill( runIndex, npro, avgDcaXYpro, 1);

      
    }

  TH2F * h1 = (TH2F*) th2f["tpc"]->Clone();
  TH2F * h2 = (TH2F*) th2f["tpc_tofmatch"]->Clone();

  h2->Divide(h1);
  h2->SetName("tof_eff");
  th2f["tof_eff"] = h2;  
  
  outfile->cd();
  
  for (auto &m : th1f)
    {
      m.second->Write();
    }

  for (auto &m : th2f)
    {
      m.second->Write();
    }

  for (auto &m : prof)
    {
      m.second->Write();
    }

  for (auto &m : prof2d)
    {
      m.second->Write();
    }


  ofstream runIndexFile;
  runIndexFile.open("RunIndexToRunNumber.txt");
  for ( auto &v : runNumVec )
    {
      runIndexFile << GetRunIndex(runNumVec,v) << ", " <<  v << "\n";
    }
  runIndexFile.close();
  
  return 0;

}

