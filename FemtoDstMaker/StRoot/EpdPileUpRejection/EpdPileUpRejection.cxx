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
#include "TH1F.h"
#include "TH2F.h"

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

#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "phys_constants.h"

#include "StRoot/EpdPileUpRejection/EpdPileUpRejection.h"

ClassImp(EpdPileUpRejection);

//__________________________________________________________________________ 
EpdPileUpRejection::EpdPileUpRejection()
{

  mEpdGeom = new StEpdGeom();  
  FillArrays();

}

//__________________________________________________________________________ 
EpdPileUpRejection::~EpdPileUpRejection()
{
}

//__________________________________________________________________________ 
void EpdPileUpRejection::FillArrays()
{

  std::vector<float> fiftyCut{118.75, 115.75, 114.75, 112.75, 109.75, 105.75, 102.25, 98.25, 95.25, 92.25};
  std::vector<float> twentyCut{97.75, 89.25, 87.75, 85.75, 83.25, 80.25, 77.25, 74.25, 71.75, 68.75};
  std::vector<float> tenCut{83.25, 72.75, 71.25, 69.75, 67.75, 65.25, 63.25, 61.75, 59.75, 57.25};

  std::vector<float> Sig1Cut{34.9136, 29.0437, 23.521, 18.4237, 14.6578, 11.4118, 9.30236, 8.00039, 6.24387,0};
  std::vector<float> Sig2Cut{32.1617, 26.6969, 21.1897, 16.2206, 12.9155, 9.79838, 8.06289, 6.94488, 5.22138,0};
  std::vector<float> Sig3Cut{29.4098, 24.3502, 18.8585, 14.0175, 11.1732, 8.18495, 6.82341, 5.88937, 4.19889,0};
  
  epdCutArray.push_back( fiftyCut );
  epdCutArray.push_back( twentyCut );
  epdCutArray.push_back( tenCut );

  tofCutArray.push_back( Sig1Cut );
  tofCutArray.push_back( Sig2Cut );
  tofCutArray.push_back( Sig3Cut );

}


Float_t EpdPileUpRejection::EpdNMip(StPicoDst * mPicoDst)
{


  _mPicoDst = mPicoDst;


   double countnmip = -999; 

   //-----------------Get EPD information------------------------------
   Int_t nepdHits = mPicoDst->numberOfEpdHits();
   StPicoEpdHit *epdHit;


   if(nepdHits >= 75){
     
     countnmip = 0; 
     

     for(Int_t iHit=0; iHit<nepdHits; iHit++){

       epdHit = mPicoDst->epdHit(iHit);

       if (!epdHit) continue;
       float mip = epdHit->nMIP();
       TVector3 StraightLine_center;
       StraightLine_center = mEpdGeom->TileCenter(epdHit->id()) - mPicoDst->event()->primaryVertex(); 
       double eta_epd_center = StraightLine_center.Eta();

       if( ( eta_epd_center>-4.0 && eta_epd_center<-2.5)  && (  mip > 0.3 && mip <= 4.0 ) ) countnmip += mip;

     }
     
   }

   return countnmip;

}


Int_t EpdPileUpRejection::Centrality(int val)
{
  //---------------------Centrality----------------------------------
  
   int gRefMult = val;
   int centrality;
   int centFull[11]={8,10,12,14,18,22,28,34,40,48,85};
   if      (gRefMult>=centFull[10]) centrality=-1;
   else if (gRefMult>=centFull[9]) centrality=0;
   else if (gRefMult>=centFull[8]) centrality=1;
   else if (gRefMult>=centFull[7]) centrality=2;
   else if (gRefMult>=centFull[6]) centrality=3;
   else if (gRefMult>=centFull[5]) centrality=4;
   else if (gRefMult>=centFull[4]) centrality=5;
   else if (gRefMult>=centFull[3]) centrality=6;
   else if (gRefMult>=centFull[2]) centrality=7;
   else if (gRefMult>=centFull[1]) centrality=8;
   else if (gRefMult>=centFull[0]) centrality=9;
   else centrality = -1;


   return centrality; 

}

std::pair<int,int> EpdPileUpRejection::FxtMultAndPiPDuCount(StPicoDst * mPicoDst)
{

  _mPicoDst = mPicoDst;

 StPicoEvent * picoEvent = mPicoDst->event();
   int countrefmult=0;     // for centrality determination
   int counttofpipdu = 0;
   TVector3 pVertex = picoEvent->primaryVertex();

  const Int_t nTrack = mPicoDst->numberOfTracks();

   // primary track loop for determine refmult ----------------------------------------------
   for (Int_t itr=0;itr<nTrack;itr++) {
      const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
      if(!ptrk)  continue;

      const Float_t pt  = ptrk->pMom().Perp(); // zero for global tracks
      const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
      const Int_t nHitsFit = ptrk->nHitsFit();
      const Int_t nHitsPoss = ptrk->nHitsMax();
      const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;
      
      //      if( eta < -1.5 || -0.5 < eta) continue;
      if( fabs(dca)>3.0 ) continue;
      if( nHitsFit < 15 ) continue;
      if( quality < 0.52 ) continue;

      int tofIndex = ptrk->bTofPidTraitsIndex();
      bool TofMatch = kFALSE;
      StPicoBTofPidTraits* tofPidTraits;

      if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex); //GNX 
      if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
      if (!TofMatch) continue;

      float beta = getTofBeta(mPicoDst, ptrk, pVertex);
      bool tofAvailable = !isnan(beta) && beta > 0;
      bool tofPion     = tofAvailable && isTofPion(ptrk, beta, pVertex);//this is hybrid pid, not always require tof
      bool tofProton   = tofAvailable && isTofProton(ptrk, beta, pVertex);
      bool tofDeuteron = tofAvailable && isTofDeuteron(ptrk, beta, pVertex);
      if (!tofAvailable) continue;

      if (pt > 1.5) tofPion = 0;
      if (pt > 3.0) tofProton = 0;
      if (pt > 3.0) tofDeuteron = 0;
      if (tofPion || tofProton || tofDeuteron) counttofpipdu++;
  
      countrefmult++;

   }

   std::pair<int,int> mult_pipdu_pair = make_pair(countrefmult,counttofpipdu);
   return mult_pipdu_pair;

}


void EpdPileUpRejection::PileUpEvent(StPicoDst * mPicoDst,Int_t &GoodEpd, Int_t &GoodTof)
{

  _mPicoDst = mPicoDst;

  GoodEpd = 0;
  GoodTof = 0;

  std::pair<int,int> refmult_countpipdu_pair = this->FxtMultAndPiPDuCount(mPicoDst);

  int countrefmult  = refmult_countpipdu_pair.first;
  int counttofpipdu = refmult_countpipdu_pair.second;

   int centrality = Centrality( countrefmult );
   if (centrality < 0) return;

   int countnmip = this->EpdNMip(mPicoDst);
        
   double cutEpd1 = epdCutArray.at(0).at(centrality);
   double cutTof1 = tofCutArray.at(0).at(centrality);

   double cutEpd2 = epdCutArray.at(1).at(centrality);
   double cutTof2 = tofCutArray.at(1).at(centrality);

   double cutEpd3 = epdCutArray.at(2).at(centrality);
   double cutTof3 = tofCutArray.at(2).at(centrality);
   
   if (countnmip < cutEpd3 ) GoodEpd=3; 
   if (countnmip < cutEpd2 ) GoodEpd=2; 
   if (countnmip < cutEpd1 ) GoodEpd=1; 
   if (countnmip <= 0 ) GoodEpd=0; 

   if (counttofpipdu > cutTof3 ) GoodTof=3; 
   if (counttofpipdu > cutTof2 ) GoodTof=2; 
   if (counttofpipdu > cutTof1 ) GoodTof=1; 

}


bool EpdPileUpRejection::isTofPion(StPicoTrack const* const trk, float beta, TVector3 const& vtx) const {
   bool tofPion = false;

   if (beta > 0) {
      double ptot = trk->gMom(vtx, _mPicoDst->event()->bField()).Mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < 0.02 ? true : false;
   }
   return tofPion;
}

bool EpdPileUpRejection::isTofProton(StPicoTrack const* const trk, float beta, TVector3 const& vtx) const {
   bool tofProton = false;

   if (beta > 0) {
      double ptot = trk->gMom(vtx, _mPicoDst->event()->bField()).Mag();
      float beta_pr = ptot / sqrt(ptot * ptot + M_PROTON * M_PROTON);
      tofProton = fabs(1 / beta - 1 / beta_pr) < 0.03 ? true : false;
   }
   return tofProton;
}

bool EpdPileUpRejection::isTofDeuteron(StPicoTrack const * const trk, float beta, TVector3 const& kfVtx) const {
   bool tofDeuteron = false;

   if (beta > 0) {
      double ptot = trk->gMom(kfVtx, _mPicoDst->event()->bField()).Mag();
      float beta_du = ptot / sqrt(ptot * ptot + M_DEUTERON * M_DEUTERON);
      tofDeuteron = fabs(1 / beta - 1 / beta_du) < 0.03 ? true : false;
   }

   return tofDeuteron;
}

float EpdPileUpRejection::getTofBeta(StPicoDst * mPicoDst, StPicoTrack const* const trk, TVector3 const& vtxTV3) const {


   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();


   if (index2tof >= 0) {

      StPicoBTofPidTraits const* const tofPid = mPicoDst->btofPidTraits(index2tof);


      if (tofPid) {


         beta = tofPid->btofBeta();
         if (beta < 1e-4) {
            TVector3 const btofHitPosTV3 = tofPid->btofHitPos();
            StThreeVectorF const vtx(vtxTV3.X(),vtxTV3.Y(),vtxTV3.Z());
            StThreeVectorF btofHitPos(btofHitPosTV3.X(),btofHitPosTV3.Y(),btofHitPosTV3.Z());



            StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
            float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();


            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();

         }
      }
   }
   return beta;
}
