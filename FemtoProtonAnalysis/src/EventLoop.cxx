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

  //  TH2F * h3 = new TH2F("h3","h3",500,-4,4,500,0,4);
  
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
      vector<vector<double>> qrs = make_all_q_s( track_eff_pair_vec, 4, 4);      
      cpc->FillProfile(fxt3,qrs);	  

    }

  cpc->FactorialCumulantsToNormalCumulants();
  
  return 0;

}


int EventLoopSystematic(
	      TChain * tc,
	      long int nentries,
	      InputParameterList & pl,
	      vector<CumulantProfileContainer*> &cpc,
	      ProtonEfficiency * eff
	      )
{


  TH2F * h[10];

  for (int i=0;i<10;i++)
    {
      h[i] = new TH2F(TString::Format("h_%i",i),"",500,0,2,500,0,2.5);
    }
  
  StFemtoEvent * event = new StFemtoEvent();  
  tc->SetBranchAddress("StFemtoEvent",&event);

  if (cpc.size() != 10)
    {
      cout << "Sys Array must be 10!" << endl;
      return 0;
    }
  
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
      if ( !epd_tof.first ) continue;
      
      vector<vector<pair<int,double>>> track_eff_pair_vec_sys;

      for (int i=0;i<10;i++)
	{
	  vector<pair<int,double>> temp;
	  track_eff_pair_vec_sys.push_back(temp);
	}
      
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

	  double trackEff = -1;
	  bool tofmatch =false;

	  if ( tofmass2 > pl.Mass2Low() && tofmass2 < pl.Mass2High()) tofmatch =true;

	  trackEff = eff->GetEff(charge,pt,pz,false);	      

	  //	  cout << "TofMass2 = " << tofmass2 << endl;
	  //	  cout << "TofMatch = " << tofmatch << endl;
	  
	  if (trackEff<0) continue;
	  if (charge<0)   continue;
	  
	  pair<int,double> track_eff_pair = make_pair(charge,trackEff);

	  /*
	  double dRap = 0.05;
	  for (int i=0;i<10;i++)
	    {
	      if ( rap > (pl.RapidityLow() - dRap*i) && rap < (pl.RapidityHigh() + dRap*i) )
		{
		  track_eff_pair_vec_sys[i].push_back(track_eff_pair);		  
		}
	    }
	  */	  

	  double dP=0.2;
	  
	  for (int i=0;i<10;i++)
	    {
	      double p = sqrt( pow(pt,2) + pow(pz,2) );
	      double pMax = pl.Mom() + dP*i;

	      //	      cout << "    p=" << p  << " pMax=" << pMax << " tofMatch=" << tofmatch << endl;
	      
	      if ( (p < pMax)  || tofmatch )
		{
		  //		  cout << "passed" << endl;
		  h[i]->Fill(rap,pt);
		  track_eff_pair_vec_sys[i].push_back(track_eff_pair);		  
		}
	    }


	}


      //Where the factorial moment calculation takes place

      for (int i=0;i<10;i++)
	{
	  vector<vector<double>> qrs = make_all_q_s( track_eff_pair_vec_sys[i], 4, 4);      
	  cpc[i]->FillProfile(fxt3,qrs);	  
	}
    }

  for (int i=0;i<10;i++)
    {
      cpc[i]->FactorialCumulantsToNormalCumulants();
    }


  TFile * of = new TFile("acceptanceCheck.root","recreate");
  of->cd();
  for (int i=0;i<10;i++)
    {
      h[i]->Write();
    }
  
  return 0;

}


int EventLoopBootstrap(
	      TChain * tc,
	      long int nentries,
	      InputParameterList & pl,
	      CumulantProfileContainer* cpc,
	      std::vector<CumulantProfileContainer*> cpc_vec,
	      ProtonEfficiency * eff,
	      TRandom3 * rand
	      )
{

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
	  cout << percent  << "%" << " iEntry =" << iEntry << endl;
	  percent++;
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
      if ( !epd_tof.first || !epd_tof.second ) continue;
      
      vector<pair<int,double>> track_eff_pair_vec;

      std::vector<double> avg_dcaXY;
      
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

	  TVector3 pMom = TVector3(trk.GetPx(),trk.GetPy(),trk.GetPz()).Unit();
	  float cos1 = pMom.Perp();
	  float dcaD = dcaX*pMom.Y()/cos1  - dcaY*pMom.X()/cos1;
	  int sign = (dcaD > 0 ) ? 1 : -1;
	  float dcaXY = sign*sqrt(dcaX*dcaX + dcaY*dcaY);
	  avg_dcaXY.push_back( dcaXY );

	  
	  //Make Track Cuts
	  if ( pt < pl.PtLow() ) continue;
	  if ( pt > pl.PtHigh() ) continue;
	  if ( rap < pl.RapidityLow() ) continue;
	  if ( rap > pl.RapidityHigh() ) continue;
	  if ( fabs(nsigpro) > pl.NSigmaProtonCut()) continue;
	  if ( fabs(nsigpi) < pl.NSigmaPionCut()) continue;
	  if ( nhitsfit <= pl.NHitsFitMin()) continue; 

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
	  
	  pair<int,double> track_eff_pair = make_pair(charge,trackEff);

	  track_eff_pair_vec.push_back(track_eff_pair);
	  
	}

      double avgDcaXY = accumulate(avg_dcaXY.begin(),avg_dcaXY.end(),0.0)/((float)avg_dcaXY.size());
      if ( avgDcaXY > 0.4 ) continue;
      
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

  cpc->FactorialCumulantsToNormalCumulants();

  for (auto &v : cpc_vec )
    {
      v->FactorialCumulantsToNormalCumulants();
    }
      
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
  
  for (int iCent=minCent;iCent<=maxCent;iCent++)
    {

      cout << "Cent=" << iCent << endl;
      
      for (int iEntry=0;iEntry<nentries_perCentBin;iEntry++)
	{     
	  vector<pair<int,double>> track_eff_pair_vec;
	  int poissonVal = rand->Poisson(iCent);
	  
	  for (int iTrack=0;iTrack<poissonVal;iTrack++)
	    {
	      if (brand->Rndm() > eff) continue;
	      pair<int,double> track_eff_pair = make_pair(1,eff);	      
	      track_eff_pair_vec.push_back(track_eff_pair);
	    }
	  
	  //Where the factorial moment calculation takes place
	  vector<vector<double>> qrs = make_all_q_s( track_eff_pair_vec, 4, 4);      
	  cpc->FillProfile(iCent,qrs);	  
	  
	}
    }

  return 0;

}
