
#include <iostream>
#include <cstdio>
#include <vector>
#include <utility>
#include "Riostream.h"

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1D.h"

#include "StMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StTriggerData.h"
#include "StTriggerIdCollection.h"
#include "StRunInfo.h"

#include "StRoot/StFemtoTrack/StFemtoTrack.h"
#include "StRoot/StFemtoEvent/StFemtoEvent.h"

#include "StRoot/StFemtoDstMaker/StFemtoDstMaker.h"

#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "phys_constants.h"

#include "StRoot/EpdPileUpRejection/EpdPileUpRejection.h"

ClassImp(StFemtoDstMaker)

StFemtoDstMaker::StFemtoDstMaker(char * name) : 
  StMaker(name){
}

StFemtoDstMaker::~StFemtoDstMaker()
{
  delete mPicoDstMaker;
} 

Int_t StFemtoDstMaker::Init()
{

  TString filename = "";
  filename.Append(mFileIndex);
  filename.Append(".fDst.root");
  filename.Prepend(mOutDir);
  mOutfile = new TFile(filename,"recreate");

  mFemtoEvent = new StFemtoEvent();
  fDstTree = new TTree("fDst","fDst");
  fDstTree->Branch("StFemtoEvent",mFemtoEvent);

  InitRunIndex();
  mEpdPileUpRejection = new EpdPileUpRejection();
  fullEventCounter=0;

  return kStOK;
}

Int_t StFemtoDstMaker::Finish() 
{
  mOutfile->cd();
  fDstTree->Write();
  return kStOK;
}

Int_t StFemtoDstMaker::Make() 
{

  mPicoDstMaker = (StPicoDstMaker *)GetMaker("PicoDst");
  if (!mPicoDstMaker){
    fputs("ERROR: StFemtoDstMaker::Init() - Can't get pointer to StPicoDstMaker!", stderr);
    return kStFATAL;
  }

  StPicoDst *mPicoDst = NULL;
  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst){
    fputs("ERROR: StFemtoDstMaker::Init() - Can't get pointer to StPicoDst!", stderr);
    return kStFATAL;
  }

  StPicoEvent * event = mPicoDst->event();

  if (!event) {
    cout << "StFemtoDstMaker::Make() No Event Found!" << endl;
    return kStOK;
  }

  Int_t runID = event->runId();
  if(! GoodRun(runID)) return kStOK;

  bool goodTrig=false;
  Int_t nTrig = static_cast<int>(event->triggerIds().size());

  for (Int_t iTrig=0; iTrig<nTrig; iTrig++)
    {
      //Check if triggerId = 620052
      if ( event->triggerIds().at(iTrig) == 620052 )
	{
	goodTrig=true;
	}

      //Check if triggerId = 620053
      if ( event->triggerIds().at(iTrig) == 620053 )
	{
	goodTrig=true;
	}

    }  
  //Skip0 all other triggers
  if (goodTrig==false) return kStOK;

  //Vertex Cut
  TVector3 pVtx = event->primaryVertex();
  Double_t vx = event->primaryVertex().X();
  Double_t vy = event->primaryVertex().Y();
  Double_t vz = event->primaryVertex().Z();
  

  //wide vertex cut
  if ( vz < 195.0 || vz > 205.0 ) return kStOK;
  double vr = TMath::Sqrt( vx*vx + (vy+2)*(vy+2) );
  if ( vr > 5 ) return kStOK;

  //Get NMIP
  Float_t epdNMip = mEpdPileUpRejection->EpdNMip(mPicoDst);

  std::pair<int,int> FxtMult_PiPDu_pair  = mEpdPileUpRejection->FxtMultAndPiPDuCount(mPicoDst);
  Int_t gMult      = FxtMult_PiPDu_pair.first;
  Int_t countpipdu = FxtMult_PiPDu_pair.second;

  mFemtoEvent->Clear();
  
  int refmult=0;
  int refmult3=0;
  double nCharge=0;

  std::vector< StFemtoTrack > trkArray;

  Double_t p_m = 0.938272; // GeV
  Double_t pi_m = 0.13957; // GeV

  Bool_t fullEvent = true;//( fullEventCounter % 10 == 0 ) ? true : false;
  fullEventCounter++;
  float mBField  = event->bField();


  for(unsigned int trackIndex=0; trackIndex< mPicoDst->numberOfTracks();trackIndex++)
    {
      
      StPicoTrack * mPicoTrack = mPicoDst->track(trackIndex);

      //Only Primary Tracks
      if ( !mPicoTrack->isPrimary() ) continue;

      //Tracks dca cut
      Float_t dca    = TMath::Abs( mPicoTrack->gDCA(vx,vy,vz) );
      Float_t dcaX  = mPicoTrack->gDCAx(vx);
      Float_t dcaY  = mPicoTrack->gDCAy(vy);
      Float_t dcaZ   = mPicoTrack->gDCAz(vz);

      //Lets check the difference
      StPicoPhysicalHelix helix_yu = mPicoTrack->helix(mBField);
      Float_t dca_yu = helix_yu.geometricSignedDistance(pVtx);
      if(dca_yu < 0) dca_yu = fabs(dca);

      //      Float_t diff = dca - dca_yu;
      //      if( fabs(diff)>0.00001) cout << "DCA = " << dca << " Yu's DCA = " << dca_yu << " DIFF = " << dca-dca_yu << endl; 


      if ( dca > 3.0 ) continue;
      
      //Remove broken tracks
      Float_t nHitsFit = mPicoTrack->nHitsFit();
      Float_t nHitsMax = mPicoTrack->nHitsMax();
      if( (nHitsFit / nHitsMax) < 0.51 ) continue;

      //Only use tracks below -3 sigma of proton band
      refmult++;

      Float_t pt = mPicoTrack->pMom().Perp();
      Float_t q = mPicoTrack->charge();
      Float_t px = mPicoTrack->pMom().X();
      Float_t py = mPicoTrack->pMom().Y();
      Float_t pz = mPicoTrack->pMom().Z();
      Float_t nSigProton = mPicoTrack->nSigmaProton();
      Float_t nSigPion = mPicoTrack->nSigmaPion();

      Float_t protonRapidity = fabs( ParticleRapidity( pt, pz, p_m ) );
      Float_t pionRapidity   = fabs( ParticleRapidity( pt, pz, pi_m ) );
      Float_t tofProtonMass =-999;
      Float_t eta = 0.5 * TMath::Log( (mPicoTrack->pMom().Mag() + mPicoTrack->pMom().Z()) / (mPicoTrack->pMom().Mag() - mPicoTrack->pMom().Z()) );


      if ( mPicoTrack->nSigmaProton() < -3 )
	{
	  refmult3++;
	  double mom = sqrt( px*px + py*py + pz*pz );
	  if ( mom > 0.2 && mom < 1.0 && eta < -0.01 && eta > -2.0)
	    {
	      double eff = GetPionEfficiency( pionRapidity, sqrt(px*px + py*py) );
	      nCharge += ( eff > 0 ) ? 1.0/eff : 0; 
	    }
	} 


      Float_t T0=-999;
      Float_t T0_p=-999;
      Float_t T0_pi=-999;



      //Check if full event is turned off and then exclude non proton particles
      if (!fullEvent && (fabs(nSigProton) > 3.2) ) continue;
      
      if (fabs(nSigProton) < 2 || fabs(nSigPion) < 2)
	{
	  if (fabs(nSigProton) < 2 ) T0_p  = this->GetT0(mPicoTrack,event->primaryVertex(),p_m);
	  if (fabs(nSigPion) < 2) T0_pi = this->GetT0(mPicoTrack,event->primaryVertex(),pi_m);
	}

      T0 = (fabs(T0_pi) < fabs(T0_p)) ? T0_pi : T0_p;

      //Fill Mass if availabe
      Int_t tofPidIndex = mPicoTrack->bTofPidTraitsIndex();
      //Check if track has a btof info
      if (tofPidIndex >= 0)
	{
	  StPicoBTofPidTraits * pidTraits = mPicoDst->btofPidTraits( tofPidIndex );
	  
	  //btof quality cuts
	  if (pidTraits->btofMatchFlag() > 0)
	    {
	      if (pidTraits->btofBeta() > 0.0 && 
		  TMath::Abs(pidTraits->btofYLocal()) < 1.6 && 
		  TMath::Abs(pidTraits->btofZLocal()) < 3.0 )
		{
		  double beta = pidTraits->btofBeta();
		  double momMag = mPicoTrack->pMom().Mag();
		  double mass2 = ( momMag*momMag*( 1 - beta*beta ) / ( beta*beta ) );

		  tofProtonMass=TMath::Sqrt(mass2);
		}
	    }
	}

      StFemtoTrack femtoTrack;

      femtoTrack.SetPx(px);
      femtoTrack.SetPy(py);
      femtoTrack.SetPz(pz);
      femtoTrack.SetNHitsFitAndCharge( q, nHitsFit);
      femtoTrack.SetDcaX(dcaX);
      femtoTrack.SetDcaY(dcaY);
      femtoTrack.SetDcaZ(dcaZ);
      femtoTrack.SetNSigmaPion(nSigPion);
      femtoTrack.SetNSigmaProton(nSigProton);
      femtoTrack.SetTofMass(tofProtonMass);
      femtoTrack.SetTofT0(T0);

      trkArray.push_back(femtoTrack);

    }

  cout << " nCharge = " << nCharge << endl;

  //Fill Event Varialbles 
  mFemtoEvent->SetRunID(runID);
  mFemtoEvent->SetFxtMult(refmult);
  mFemtoEvent->SetFxtMultTofMatch(nCharge);
  mFemtoEvent->SetFxtMult3(refmult3);
  
  mFemtoEvent->SetVr(vr);
  mFemtoEvent->SetVz(vz);
  
  mFemtoEvent->SetEpdNMip(epdNMip);
  mFemtoEvent->SetPiPDu(countpipdu);
  
  mFemtoEvent->SetStFemtoTrackArray(trkArray);
  //  mFemtoEvent->PrintEvent();
  
  //  mFemtoEvent->SetFullEvent(fullEvent);

  fDstTree->Fill();

  return kStOK;  

}

void StFemtoDstMaker::InitRunIndex(){

  ifstream in;
  in.open(runNumFile.Data());
  
  while(1){
    int runNum = -999;
    in >> runNum;
    if (runNum > 0) cout << "StMomentQAMaker::InitRunIndex() Adding run index: "<< runNum << endl;
    if (runNum > 0) run_vec.push_back(runNum);
    if (!in.good()) break;
  }

  std::sort( run_vec.begin(),run_vec.end() );
  
}

int StFemtoDstMaker::GetRunIndex(int runId){

  int index = -999;
  
  for(std::vector<int>::iterator iRun = run_vec.begin() ; iRun != run_vec.end(); iRun++)
    {

    if ( runId == (*iRun) ) {
      index = std::distance(run_vec.begin(),iRun); 
      break;
    }    

  }

  if ( index >= 0) return index;
  else {
    cout << "bad run..." << endl;
    return index;
  }
}

bool StFemtoDstMaker::GoodRun(int runId)
{
  return ( GetRunIndex(runId) >= 0 ) ? true : false;
}


Float_t StFemtoDstMaker::GetT0(StPicoTrack const* const trk, TVector3 const& vtxTV3,float mass) const {
   
  int index2tof = trk->bTofPidTraitsIndex();

   float beta1 = std::numeric_limits<float>::quiet_NaN();
   float beta2 = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0) {
      StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if (tofPid) {
	beta1 = tofPid->btof();
	TVector3 const btofHitPosTV3 = tofPid->btofHitPos();
	StThreeVectorF const vtx(vtxTV3.X(),vtxTV3.Y(),vtxTV3.Z());
	StThreeVectorF btofHitPos(btofHitPosTV3.X(),btofHitPosTV3.Y(),btofHitPosTV3.Z());
	
	StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
	float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
	//	float tof = tofPid->btof();
	float alpha = trk->pMom().Mag()/mass;
	float v = (C_C_LIGHT/1.e9) * TMath::Sqrt( alpha*alpha / ( 1 + alpha*alpha ) );
	beta2 = L/v;
      }
   }
   
   //   cout << "beta1=" << beta1 << endl;
   //   cout << "beta2=" << beta2 << endl;

   return beta1 - beta2;
}

Double_t StFemtoDstMaker::ParticleRapidity( Double_t pt, Double_t pz, Double_t mass )
{
  Double_t energy = sqrt( pow(pt,2) + pow(pz,2) + pow(mass,2) );
  return 0.5 * log( ( energy + pz) / (energy - pz) );
}


int StFemtoDstMaker::LoadTPCEff(TString filename)
{

  TFile * f = new TFile(filename,"read");

  std::vector<TH1D*> h_vec(100,(TH1D*) NULL);
  std::vector<TF1*> f_vec(100,(TF1*) NULL);

  cout << "LoadTPCEff() loading histograms and fits ..." << endl;
  for (int i=3;i<97;i++) 
    {
      h_vec[i] = (TH1D*) f->Get(TString::Format("slice_bin%i",i).Data());
      f_vec[i] = (TF1*) f->Get(TString::Format("f_slice_bin%i",i).Data());
    }

  cout << " hists and fits loaded" << endl;

  _h1d_eff = h_vec;
  _tf1_eff = f_vec;

  return 0;

}

int StFemtoDstMaker::GetPionEffBin( double y )
{
  float delta=0.02;
  float start=3*0.02;

  for ( int i =3; i < 97; i++)
    {
      if ( y >= start && y < start+delta ) return i;
      else start += delta;
    }
  
  return -1;

}

double StFemtoDstMaker::GetPionEfficiency(double y, double pt)
{

  //Get the index for the correct histogram and fit
  int index = GetPionEffBin( y );
  if (index<0) return -1;

  //  cout << "index=" << index << endl;

  TH1D * h = _h1d_eff[index];
  TF1 * f = _tf1_eff[index];
  
  if ( !h ) cout << "hist doesnt exist" << endl;
  if ( !f ) cout << "fit doesnt exist" << endl;

  float xMin = h->GetBinLowEdge(h->FindFirstBinAbove(0.01));
  float xMax = 2.0;
  if ( pt <= xMin || pt >= xMax ) return -1;

  double eff = f->Eval( pt );

  if ( eff > 1.0 || eff < 0.01 ) return -1; //Set to 1% cut off
  
  return eff;
}


