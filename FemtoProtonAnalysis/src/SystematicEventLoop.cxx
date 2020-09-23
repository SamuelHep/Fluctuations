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
#include "SystematicEventLoop.h"

#include "analysisUtil.h"

using namespace std;


/*************
This is similar to the event loop but does not calcuate the cumulants, only the raw moments
**************/

void GetProfilesAndCalculateSystematics(TString infilename,TString outfilename, int nBootstraps, TString sysLabel, double start, double sysDelta)
{

  TFile * f = new TFile(infilename,"read");
  std::vector< CumulantProfileContainer* > cpc_vec;  
  std::vector<TGraphErrors*> compare_graphs;

  //Get the profiles from root file
  int iContainer=0;
  for (int iSys=0;iSys<10;iSys++)
    {
      for ( int iBootstrap=0;iBootstrap<nBootstraps;iBootstrap++)
	{
	  cpc_vec.push_back(new CumulantProfileContainer(f,iContainer));
	  iContainer++;
	}
    }


  //Bin labels to apply
  std::vector<double> binLabels = {118,170,240,308 ,395};
  std::vector<int> binEdges =  {12,19,27,39,47 ,80};
  
  //Calculate the factorial moments and rebin
  for ( auto &c : cpc_vec )
    {
      c->FactorialCumulantsToNormalCumulants();
      c->ReBinAllGraphs(binEdges,binLabels);
    }


  //Calculate the bootstrap errors for each systematic
  std::vector< CumulantProfileContainer* > cpc_vec_primary;

  iContainer=0;
  for (int iSys=0;iSys<10;iSys++)
    {
      int index_of_primary = 0;
      std::vector< CumulantProfileContainer* > cpc_vec_bs;
      
      for ( int iBootstrap=0;iBootstrap<nBootstraps;iBootstrap++)
	{
	  if ( iBootstrap==0)
	    {
	      index_of_primary=iContainer;
	    }
	  else
	    {
	      cpc_vec_bs.push_back( cpc_vec[iContainer] );
	    }
	  iContainer++;
	}
      ComputeBootstrapErrors(cpc_vec[index_of_primary],cpc_vec_bs, compare_graphs);
      cpc_vec_primary.push_back(cpc_vec[index_of_primary]);
    }


  vector<TGraphErrors*> comp_plots = MakeSysComparisonPlots(
							    cpc_vec_primary[0],
							    cpc_vec_primary,
							    sysLabel,
							    start,
							    sysDelta
							    );



  //Save all of the graphs
  TFile * fout = new TFile(outfilename,"recreate");
  fout->cd();

  cout << "Compare Graphs # = " << compare_graphs.size() << endl;

  for (auto &v : compare_graphs)
    {
      cout << "here writing" << endl;
      v->Write();
    }
  
  for (auto &v : cpc_vec_primary)
    {
      std::vector<TGraphErrors*> sys_ngr_vec = v->GetGraphVec();
      for (auto &tg : sys_ngr_vec)
	{
	  tg->Write();
	}
    }

  for (auto &v : comp_plots)
    {
      v->Write();
    }
  

}


int SystematicEventLoopPrintToFile(
				   TChain * tc,
				   long int nentries,
				   InputParameterList & pl,
				   ProtonEfficiency * eff,
				   TRandom3 * rand,
				   int nBootstraps,
				   TString sysLabel,
				   TString outfileName
			       )
{
  
  //Make CumulantContainer Object to store result
  std::vector< CumulantProfileContainer* > cpc_vec;

  int profIndex=0;

  for (int iSys=0;iSys<10;iSys++)
    {
      for (int iBootstraps=0;iBootstraps<nBootstraps;iBootstraps++)
	{
	  cpc_vec.push_back(new CumulantProfileContainer(profIndex));
      profIndex++;
	}
    }

  SysEventLoop( tc , nentries , pl , cpc_vec , nBootstraps, sysLabel,eff , rand); //Only Calculates the moment profiles
  TFile * outfile = new TFile(outfileName,"recreate");
  outfile->cd();

  profIndex=0;

  // Save bootstrap profiles
  for (int iSys=0;iSys<10;iSys++)
    {
      for ( int iBootstraps=0;iBootstraps<nBootstraps;iBootstraps++)
	{
	  vector<TProfile*> prof_vec_bs = cpc_vec[profIndex]->GetProfileVec();
	  for (auto &p : prof_vec_bs)
	    {
	      p->Write();
	    }
	  profIndex++;
	}

    }

  outfile->Close();

}



int SysEventLoop(
		 TChain * tc,
		 long int nentries,
		 InputParameterList & pl,
		 std::vector<CumulantProfileContainer*> cpc_vec,
		 int nBootstraps,
		 TString sysLabel,
		 ProtonEfficiency * eff,
		 TRandom3 * rand
		 )
{

  double pt_min_list[10] = {0.4,   0.545, 0.715, 0.915, 0.4, 0.8, 1.2, 1.6, 0.4,   0.705};
  double pt_max_list[10] = {0.545, 0.705, 0.915, 2.0  , 0.8, 1.2, 1.6, 2.0, 0.705, 2.0};

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
      double fxttof = event->GetFxtMultTofMatch();      

      //Make Event Cuts
      if ( vz > pl.VzMax() ) continue;
      if ( vz < pl.VzMin() ) continue;
      if ( vr > pl.VrMax() ) continue;
      double nmip = event->GetNMip();
      double pipdu = event->GetPiPDu();
      pair<bool,bool> epd_tof = PileUpBadEvent(fxttof,(int)fxt,nmip,pipdu);
      
      if ( sysLabel == "pileup" )
	{}
      else
	{
	  //	  if ( !epd_tof.first || !epd_tof.second ) continue;
	}

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
	  if ( sysLabel == "pt" ){}
	  else
	    {
	    if ( pt < pl.PtLow() ) continue;
	    if ( pt > pl.PtHigh() ) continue;
	    }
	  if ( sysLabel == "rap" || sysLabel == "rapdelta" ){}
	  else
	    {
	      if ( rap < pl.RapidityLow() ) continue;
	      if ( rap > pl.RapidityHigh() ) continue;
	    }
	  if ( fabs(nsigpro) > pl.NSigmaProtonCut()) continue;
	  if ( fabs(nsigpi) < pl.NSigmaPionCut()) continue;
	  if ( nhitsfit <= pl.NHitsFitMin()) continue; 

	  //if pass all cuts, make a pair of the charge and the effficiency
	  double trackEff = -1;
	  bool tofmatch =false;
	  if ( tofmass2 > pl.Mass2Low() && tofmass2 < pl.Mass2High()) tofmatch =true;

	  if ( sysLabel == "mom" )
	    {}
	  else
	    {
	      if ( sqrt(pt*pt + pz*pz) > pl.Mom() && tofmatch==false ) continue;
	    }

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
	  
	  if ( sysLabel == "pileup" )
	    {
	      for (int i=0;i<10;i++)
		{
		  track_eff_pair_vec_sys[i].push_back(track_eff_pair);
		}
	    }
	  

	  if ( sysLabel == "pt" ) 
	    {
	      for (int i=0;i<10;i++)
		{
		  
		  double pt_min = pt_min_list[i];
		  double pt_max = pt_max_list[i];
		  
		  if ( pt > pt_min && pt < pt_max )
		    {
		      track_eff_pair_vec_sys[i].push_back(track_eff_pair);		  
		    }
		}	      
	      
	    }
	  
	  if ( sysLabel == "mom" )
	    {
	      double dP=0.2;
	      
	      for (int i=0;i<10;i++)
		{
		  double p = sqrt( pow(pt,2) + pow(pz,2) );
		  double pMax = 0.5 + dP*i;
		  
		  if ( (p < pMax)  || tofmatch )
		    {
		      track_eff_pair_vec_sys[i].push_back(track_eff_pair);		  
		    }
		}	      
	    }
	  
	  if ( sysLabel == "rap" )
	    {
	      double dRap=0.05;
	      
	      for (int i=0;i<10;i++)
		{
		  //double rapMax = 0.1 + dRap*i;
		  double rapMax = 0.1;
		  double rapMin = -0.1 - dRap*i;
		  
		  if ( ((rap -1.04) < rapMax) && ((rap -1.04) > rapMin ) )
		    {
		      track_eff_pair_vec_sys[i].push_back(track_eff_pair);		  
		    }
		}	      
	    }
	  
	  if ( sysLabel == "rapdelta" )
	    {
	      double dRap=0.05;
	      
	      for (int i=0;i<10;i++)
		{
		  double rapMax = pl.RapidityHigh() - dRap*i;
		  double rapMin = pl.RapidityLow() - dRap*i;
		  
		  if ( rap < rapMax && rap > rapMin )
		    {
		      track_eff_pair_vec_sys[i].push_back(track_eff_pair);		  
		    }
		}	      
	    }	 
	}

      //End Track Loop

    
      int profIndex=0;

      for (int i=0;i<10;i++)
	{

	  //If the epd/tof pileup cut, do the cut for each sys index
	  bool skip = false;
	  if ( sysLabel == "pileup" )
	    {
	      pair<bool,bool> sys_epd_tof = PileUpBadEventVariable((int)fxt,nmip,pipdu,i);
	      if ( !sys_epd_tof.first || !sys_epd_tof.second ) skip = true;
	    }
	  
	  for ( int iBoot=0; iBoot<nBootstraps;iBoot++)
	    {
	      if ( !skip )
		{
		  vector<vector<double>> qrs = make_all_q_s( track_eff_pair_vec_sys[i], 4, 4);      
		  if (iBoot==0) cpc_vec[profIndex]->FillProfile(fxt3,qrs);	  
		  else 
		    {
		      int nFill = rand->Poisson(1);
		      for (int iFill=0;iFill<nFill;iFill++)
			{
			  cpc_vec[profIndex]->FillProfile(fxt3,qrs);
			}
		    }
		}
	      profIndex++;
	    }
	}  

    }
  
  return 0;

}


