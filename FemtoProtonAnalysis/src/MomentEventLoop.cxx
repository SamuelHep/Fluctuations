#include <iostream>
#include <vector>
#include <utility>
#include <numeric>

#include "TFile.h"
#include "TChain.h"
#include "TRandom3.h"

#include "StFemtoTrack.h"
#include "StFemtoEvent.h"

#include "InputParameterList.h"
#include "CumulantProfileContainer.h"
#include "MomentFunctions.h"
#include "ProtonEfficiency.h"
#include "MomentEventLoop.h"
#include "analysisUtil.h"

using namespace std;


/*************
This is similar to the event loop but does not calcuate the cumulants, only the raw moments
**************/
int MomentEventLoopLocal(
			 TChain * tc,
			 long int nentries,
			 InputParameterList & pl,
			 ProtonEfficiency * eff,
			 TRandom3 * rand,
			 int nBootstraps,
			 TString outfileName
			 )
{
 
  //Bin Labels and edges
  std::vector<double> binLabels = { 10, 21, 32, 50, 74, 105, 149, 251}; // <Npart>
  std::vector<int> binEdges     = {2, 6, 8, 14, 20, 30, 42, 52, 90};    //FxtMult3 Edges

  //Make CumulantContainer Object to store result
  CumulantProfileContainer * cpc = new CumulantProfileContainer();

  //Used for bootstraps
  std::vector< CumulantProfileContainer* > cpc_vec;
  for (int iBootstraps=0;iBootstraps<nBootstraps;iBootstraps++)
    {
      cpc_vec.push_back(new CumulantProfileContainer(iBootstraps));
    }
  
  MomentEventLoop( tc , nentries , pl , cpc, cpc_vec , eff , rand); //Only Calculates the moment profiles

  //Rebin the graphs... Also calculate cumulants
  cpc->ReBinAllGraphs(binEdges,binLabels);  

  //Do the same thing for the bootstraps
  for( int iBs=0;iBs<nBootstraps;iBs++ )
    {
      cpc_vec[iBs]->ReBinAllGraphs(binEdges,binLabels);  
    }

  //Compute bootstraps
  std::vector<TGraphErrors*> compare_graphs;
  ComputeBootstrapErrors(cpc,cpc_vec, compare_graphs);

  TFile * outfile = new TFile(outfileName,"recreate");
  outfile->cd();

  //Save each graph
  for( int i=0;i<cpc->NGraphs();i++)
    {
      TGraphErrors * gr = cpc->GetNGraph(i);
      gr->Write();
      delete gr;
    }

  outfile->Close();

}

int MomentEventLoopPrintToFile(
			       TChain * tc,
			       long int nentries,
			       InputParameterList & pl,
			       ProtonEfficiency * eff,
			       TRandom3 * rand,
			       int nBootstraps,
			       TString outfileName
			       )
{
  
  //Make CumulantContainer Object to store result
  CumulantProfileContainer * cpc = new CumulantProfileContainer();
  std::vector< CumulantProfileContainer* > cpc_vec;

  for (int iBootstraps=0;iBootstraps<nBootstraps;iBootstraps++)
    {
      cpc_vec.push_back(new CumulantProfileContainer(iBootstraps));
    }
  
  MomentEventLoop( tc , nentries , pl , cpc, cpc_vec , eff , rand); //Only Calculates the moment profiles

  TFile * outfile = new TFile(outfileName,"recreate");
  outfile->cd();

  // Save primary profile vec
  vector<TProfile*> prof_vec = cpc->GetProfileVec();
  for (auto &p : prof_vec)
    {
      p->Write();
    }

  // Save bootstrap profiles
  for ( int iBootstraps=0;iBootstraps<nBootstraps;iBootstraps++)
    {
      vector<TProfile*> prof_vec_bs = cpc_vec[iBootstraps]->GetProfileVec();
      for (auto &p : prof_vec_bs)
	{
	  p->Write();
	}
    }

  outfile->Close();

}


int MomentEventLoop(
		    TChain * tc,
		    long int nentries,
		    InputParameterList & pl,
		    CumulantProfileContainer* cpc,
		    std::vector<CumulantProfileContainer*> cpc_vec,
		    ProtonEfficiency * eff,
		    TRandom3 * rand
		    )
{

  //Checks
  TH1D * h_protons[300];
  for ( int i=0;i<300;i++)
    {
      h_protons[i] = new TH1D(TString::Format("h_protons_FxtMult3_%i",i),"",300,-0.5,299.5);
    }
  TH1D * h_mult = new TH1D("h_mult","",300,-0.5,299.5);
  TH1D * h_p = new TH1D("h_p","",300,-0.5,299.5);

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
      
      //Make Event Cuts
      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
      //    pair<bool,bool> epd_tof = PileUpBadEvent((int)fxt,nmip,pipdu);
      //    if ( !epd_tof.first || !epd_tof.second ) continue;

      //      if (fxt >= 210) continue;
     
      //Track Vector to be filled 
      vector<pair<int,double>> track_eff_pair_vec;
      std::vector<double> avg_dcaXY; //Dca Check // Not Used
      
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
	  
	  //Make Track Cuts
	  if ( pt < pl.PtLow() ) continue;
	  if ( pt > pl.PtHigh() ) continue;
	  if ( rap < pl.RapidityLow() ) continue;
	  if ( rap > pl.RapidityHigh() ) continue;
	  if ( fabs(nsigpro) > pl.NSigmaProtonCut()) continue;
	  if ( fabs(nsigpi) < pl.NSigmaPionCut()) continue;
	  if ( nhitsfit < pl.NHitsFitMin()) continue; 
	  if ( dca > pl.Dca() ) continue;

	  //if pass all cuts, make a pair of the charge and the effficiency
	  double trackEff = -1;
	  bool tofmatch =false;
	  if ( tofmass2 > pl.Mass2Low() && tofmass2 < pl.Mass2High()) tofmatch =true;

	  if ( sqrt(pt*pt + pz*pz) > pl.Mom() && tofmatch==false ) continue;

	  if ( sqrt(pt*pt + pz*pz) > pl.Mom() )
	    {
	      trackEff = eff->GetEff(charge,pt,pz,true);
	    }
	  else
	    {
	      trackEff = eff->GetEff(charge,pt,pz,false);
	    }

	  if (trackEff<0) continue;
	  if (charge<0)   continue;
	

	  trackEff = trackEff*pl.EffMultiplier(); //use for systematic check

	  //Holds track charge and efficiency. Used in q(r,s)
	  pair<int,double> track_eff_pair = make_pair(charge,trackEff);
	  track_eff_pair_vec.push_back(track_eff_pair);
	  
	}
    
      /*
      //Check

      if ( fxt3 < 300 )
	{
	  int Np = track_eff_pair_vec.size();
	  h_mult->Fill( (int) fxt3 );
	  h_p->SetBinContent( h_p->FindBin( fxt3 ) , h_p->GetBinContent(  h_p->FindBin( fxt3 ) ) + Np );
	  h_protons[ (int) fxt3 ]->Fill( (int) track_eff_pair_vec.size() );
	}
      */

      //Where the factorial moment calculation takes place
      vector<vector<double>> qrs = make_all_q_s( track_eff_pair_vec, 4, 4);
      
      //For boostrap, the probability that the event is filled N times is poisson.
      //Draw from poission distribution and fill N times
      int nFill = 0;

      cpc->FillProfile(fxt3,qrs);
      
      for( auto &c : cpc_vec )
	{
	  nFill = rand->Poisson(1);
	  for (int iFill=0;iFill<nFill;iFill++)
	    {
	      c->FillProfile(fxt3,qrs);
	    }
	}
    }

  /*      
  //Check
  for ( int i=0;i<300;i++)
    {
      if ( h_mult->GetBinContent(i) > 0 ) h_p->SetBinContent( i ,1.0*h_p->GetBinContent(i) / h_mult->GetBinContent(i) );
    }

  TFile * proton_hist_file = new TFile("proton_hists.root","recreate");
  for ( int i=0;i<300;i++)
    {
      if (h_protons[i]->GetEntries() >0 ) h_protons[i]->Write();
    }
  
  h_mult->Write();
  h_p->Write();

  proton_hist_file->Close();
  */

  return 0;

}


int PileUpEventLoop(
		    TChain * tc,
		    long int nentries,
		    InputParameterList & pl,
		    TString outfilename
		    )
{

  map<TString,TH2F*> th2f;

  th2f["refmult_nmip"] = new TH2F("refmult_mnmip",";FxtMult;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["gmult_nmip"] = new TH2F("gmult_mnmip",";FxtMult TofMatch;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["refmult3_nmip"] = new TH2F("refmult3_mnmip",";FxtMult3;EPD #Sigma nMip",300,0,300,200,0,200);

  th2f["refmult_nmip_cut"] = new TH2F("refmult_mnmip_cut",";FxtMult;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["gmult_nmip_cut"] = new TH2F("gmult_mnmip_cut",";FxtMult TofMatch;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["refmult3_nmip_cut"] = new TH2F("refmult3_mnmip_cut",";FxtMult3;EPD #Sigma nMip",300,0,300,200,0,200);

  th2f["refmult_pipdu"] = new TH2F("refmult_pipdu",";FxtMult;TOF #pi + p + du",300,0,300,300,0,300);
  th2f["refmult3_pipdu"] = new TH2F("refmult3_pipdu",";FxtMult3;TOF #pi + p + du",300,0,300,300,0,300);

  th2f["refmult_pipdu_cut"] = new TH2F("refmult_pipdu_cut",";FxtMult;TOF #pi + p + du",300,0,300,300,0,300);
  th2f["refmult3_pipdu_cut"] = new TH2F("refmult3_pipdu_cut",";FxtMult3;TOF #pi + p + du",300,0,300,300,0,300);

  th2f["refmult_pipdu_ebad"] = new TH2F("refmult_pipdu_ebad",";FxtMult;TOF #pi + p + du",300,0,300,300,0,300);
  th2f["refmult_pipdu_egood"] = new TH2F("refmult_pipdu_egood",";FxtMult;TOF #pi + p + du",300,0,300,300,0,300);

  th2f["refmult_nmip_pbad"] = new TH2F("refmult_mnmip_pbad",";FxtMult;EPD #Sigma nMip",300,0,300,200,0,200);
  th2f["refmult_nmip_pgood"] = new TH2F("refmult_mnmip_pgood",";FxtMult;EPD #Sigma nMip",300,0,300,200,0,200);

  th2f["gmult_nmip_pbad"] = new TH2F("gmult_mnmip_pbad",";FxtMult TofMatch;EPD #Sigma nMip",300,0,300,200,0,200);  
  th2f["gmult_nmip_pgood"] = new TH2F("gmult_mnmip_pgood",";FxtMult TofMatch;EPD #Sigma nMip",300,0,300,200,0,200);  

  th2f["pileup_summary"] = new TH2F("pileup_summary","pileup_summary",10,-0.5,9.5,7,-0.5,6.5);

  StFemtoEvent * event = new StFemtoEvent();
  
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  //Proton Mass
  double mass=0.938272;

  long int onePercent = nentries/100;
  int percent = 0;
  
  for (int iEntry=0;iEntry<nentries;iEntry++)
    {

      if ( iEntry % onePercent == 0 )
	{
	  //	  cout << percent  << "%" << " iEntry =" << iEntry << endl;
	  percent++;
	}
      
      tc->GetEntry(iEntry);
      
      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      double fxttof = event->GetFxtMultTofMatch();
      
      //Make Event Cuts
      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();

      pair<bool,bool> epd_tof = PileUpBadEvent(fxttof,(int)fxt,nmip,pipdu);
      //      if ( !epd_tof.first || !epd_tof.second ) continue;

      int cent = CentBin3(fxt3);

      th2f["refmult_nmip"]->Fill( fxt, nmip);
      th2f["gmult_nmip"]->Fill( fxttof, nmip);
      th2f["refmult3_nmip"]->Fill( fxt3, nmip);
      
      th2f["refmult_pipdu"]->Fill( fxt, pipdu);
      th2f["refmult3_pipdu"]->Fill( fxt3, pipdu);

      if ( epd_tof.first )
	{
	  th2f["refmult_nmip_cut"]->Fill( fxt, nmip);
	  th2f["gmult_nmip_cut"]->Fill( fxttof, nmip);
	  th2f["refmult3_nmip_cut"]->Fill( fxt3, nmip);
	}
	 
      if ( epd_tof.second )
	{
	  th2f["refmult_pipdu_cut"]->Fill( fxt, pipdu);
	  th2f["refmult3_pipdu_cut"]->Fill( fxt3, pipdu);
	}
      
      if ( !epd_tof.first ) th2f["refmult_pipdu_ebad"]->Fill( fxt, pipdu);
      if ( epd_tof.first )  th2f["refmult_pipdu_egood"]->Fill( fxt, pipdu);
      
      if ( !epd_tof.second ) th2f["refmult_nmip_pbad"]->Fill( fxt, nmip);
      if ( epd_tof.second )  th2f["refmult_nmip_pgood"]->Fill( fxt, nmip);
      
      if ( !epd_tof.second ) th2f["gmult_nmip_pbad"]->Fill( fxttof, nmip);
      if ( epd_tof.second ) th2f["gmult_nmip_pgood"]->Fill( fxttof, nmip);

      //      TString binLabels[] = {"AllEvents","EventsAfterCut","TotalCut","EPDCut","TOFCut","PiPDu==0", "OverlappingCut"};      

      if ( cent < 0 || cent > 9 ) continue;
      th2f["pileup_summary"]->Fill( cent, 0); // All Events
      if ( epd_tof.second && epd_tof.first ) th2f["pileup_summary"]->Fill( cent, 1); // EPD && TOF Cut Events
      if ( !epd_tof.second || !epd_tof.first ) th2f["pileup_summary"]->Fill( cent, 2); // EPD && TOF Cut Events
      if ( !epd_tof.first ) th2f["pileup_summary"]->Fill( cent, 3); // EPD Cut Events
      if ( !epd_tof.second ) th2f["pileup_summary"]->Fill( cent, 4); // TOF Cut Events
      if ( pipdu == 0 ) th2f["pileup_summary"]->Fill( cent, 5); // TOF Events PiPDu==0
      if ( !epd_tof.second && !epd_tof.first ) th2f["pileup_summary"]->Fill( cent, 6); // EPD && TOF Cut Events



      
    }
      

  
  TFile * outFile = new TFile(outfilename,"RECREATE");
  outFile->cd();
  for (auto &h : th2f )
    {
      h.second->Write();
    }
  outFile->Close();

  return 0;

}


