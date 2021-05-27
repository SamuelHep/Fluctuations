#include <iostream>
#include <vector>
#include <utility>
#include <numeric>
#include <random>

#include "TFile.h"
#include "TChain.h"
#include "TRandom3.h"

#include "StFemtoTrack.h"
#include "StFemtoEvent.h"

#include "InputParameterList.h"
#include "CumulantProfileContainer.h"
#include "MomentFunctions.h"
#include "ProtonEfficiency.h"

#include "analysisUtil.h"

using namespace std;

int EventLoop(
	      TChain * tc,
	      long int nentries,
	      InputParameterList & pl,
	      CumulantProfileContainer * cpc,
	      ProtonEfficiency * eff
	      )
{
  
  StFemtoEvent * event = new StFemtoEvent();  
  tc->SetBranchAddress("StFemtoEvent",&event);
  
  //Proton Mass
  double mass=0.938272;

  long int onePercent = nentries/10;
  int percent = 0;
  
  cout << "Total Number of entries =" << nentries << endl; 
  
  for (int iEntry=0;iEntry<nentries;iEntry++)
    {

      if ( iEntry % onePercent == 0  )
	{
	  cout << percent  << "%" << " iEntry =" << iEntry << endl;
	  percent = percent + 10;
	}
      
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
      pair<bool,bool> epd_tof = PileUpBadEvent((int)fxt,nmip,pipdu);
      //      if ( !epd_tof.first ) continue;
      
      vector<pair<int,double>> track_eff_pair_vec;
      
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

	  double nsigpro = trk.GetNSigmaProton();
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
	  if ( nhitsfit <= pl.NHitsFitMin()) continue; 

	  //	  h3->Fill(rap,pt);

	  double trackEff = -1;
	  bool tofmatch =false;
	  if ( tofmass2 > pl.Mass2Low() && tofmass2 < pl.Mass2High()) tofmatch =true;

	  if ( sqrt(pt*pt + pz*pz) > pl.Mom() &&tofmatch==false ) continue;

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

	  pair<int,double> track_eff_pair = make_pair(charge,trackEff);
	  //push the pair to the track_eff array

	  track_eff_pair_vec.push_back(track_eff_pair);
	  
	}

      //Where the factorial moment calculation takes place
      vector<vector<long double>> qrs = make_all_q_s( track_eff_pair_vec, 6, 6);      
      cpc->FillProfile(fxt3,qrs);	  

    }

  cpc->MomentsToCumulants();
  
  return 0;

}


int EventLoopSimPoisson(
	      long int nentries_perCentBin,
	      CumulantProfileContainer * cpc,
	      TRandom3 * rand,
	      int minCent,
	      int maxCent,
	      double eff
			)
{

  //Proton Mass
  double mass=0.938272;  

  TRandom3 * brand= new TRandom3();

  std::default_random_engine generator;
  std::poisson_distribution<int> distribution(1.0);

  
  for (int iCent=minCent;iCent<=maxCent;iCent++)
    {

      cout << "Cent=" << iCent << endl;
      
      for (int iEntry=0;iEntry<nentries_perCentBin;iEntry++)
	{     
	  vector<pair<int,double>> track_eff_pair_vec;
	  //int poissonVal = rand->Poisson(iCent);
	  int poissonVal =  distribution(generator);
	  	  
	  for (int iTrack=0;iTrack<poissonVal;iTrack++)
	    {
	      if (brand->Rndm() > eff) continue;
	      pair<int,double> track_eff_pair = make_pair(1,eff);	      
	      track_eff_pair_vec.push_back(track_eff_pair);
	    }
	  
	  //Where the factorial moment calculation takes place
	  vector<vector<long double>> qrs = make_all_q_s( track_eff_pair_vec, 6, 6);      
	  cpc->FillProfile(iCent,qrs);	  
	  
	}
    }

  cpc->MomentsToCumulants();

  return 0;

}


void RunSim(TString filename)
{

  CumulantProfileContainer * cpc = new CumulantProfileContainer();
  TRandom3 * rand = new TRandom3(123);
  
  std::default_random_engine generator;
  std::poisson_distribution<int> distribution(4.1);

  EventLoopSimPoisson( 100000,
		       cpc,
		       rand, // Random
		       1,2, // Min, Max
		       0.5
		       );

  TFile * outfile = new TFile( filename, "recreate");
  outfile->cd();

  for ( auto & key_graph : cpc->GetGraphMap() )
    {
      key_graph.second->Write();
    }

  for ( auto & key_profile : cpc->GetProfileMap() )
    {
      key_profile.second->Write();
    }

  outfile->Close();

}
