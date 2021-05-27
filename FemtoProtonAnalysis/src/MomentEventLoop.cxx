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

  cpc->MomentsToCumulants();

  //Rebin the graphs... Also calculate cumulants
  cpc->ReBinAllGraphs(binEdges,binLabels);  

  //Do the same thing for the bootstraps
  for( int iBs=0;iBs<nBootstraps;iBs++ )
    {
      cpc_vec[iBs]->MomentsToCumulants();
      cpc_vec[iBs]->ReBinAllGraphs(binEdges,binLabels);  
    }

  //Compute bootstraps
  ComputeBootstrapErrors(cpc,cpc_vec);

  TFile * outfile = new TFile(outfileName,"recreate");
  outfile->cd();

  //Save each graph  
  for ( auto & key_graph : cpc->GetGraphMap() )
    {
      key_graph.second->Write();
    }


  for( int iBs=0;iBs<nBootstraps;iBs++ )
    {
      for ( auto & key_graph : cpc_vec[iBs]->GetGraphMap() )
	{
	  key_graph.second->Write();
	}
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

  // Save primary profiles
  for (auto &key_profile : cpc->GetProfileMap() )
    {
      key_profile.second->Write();
    }

  // Save bootstrap profiles
  for ( int iBootstraps=0;iBootstraps<nBootstraps;iBootstraps++)
    {
      for (auto &key_profile : cpc_vec[iBootstraps]->GetProfileMap() )
	{
	  key_profile.second->Write();
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
      
      //      cout << iEntry << endl;
      tc->GetEntry(iEntry);
      
      //Get Event Variables
      double vz   = event->GetVz();
      double vr   = event->GetVxy();
      double fxt3 = event->GetFxtMult3();
      double fxt  = event->GetFxtMult();
      double fxttof = event->GetFxtMultTofMatch();      
      
      //      double centBin3 = CentBin3( fxt3 );
      //      if ( centBin3 < 0 ) continue;
      //Make Event Cuts
      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;

      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
     
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
    
      //Where the factorial moment calculation takes place
      vector<vector<long double>> qrs = make_all_q_s( track_eff_pair_vec, 6, 6);
      
      //For boostrap, the probability that the event is filled N times is poisson.
      //Draw from poission distribution and fill N times
      int nFill = 0;

      cpc->FillProfile( (int)fxt3,qrs);
      
      for( auto &c : cpc_vec )
	{
	  nFill = rand->Poisson(1);
	  for (int iFill=0;iFill<nFill;iFill++)
	    {
	      c->FillProfile( (int) fxt3,qrs);
	    }
	}
    }

  return 0;

}





