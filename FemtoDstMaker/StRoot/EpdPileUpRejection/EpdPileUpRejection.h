#ifndef STAR_REJECTION_MAKER
#define STAR_REJECTION_MAKER

#ifndef StMaker_H
#include "StMaker.h"
#include "TVector3.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#endif

#include <map>
#include <utility>

//Classes Used in the .cxx File                                              
class TFile;
class TBranch;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class TH2F;
class TH1F;
class TofPid;


class EpdPileUpRejection : public TObject {
  
 public:
  EpdPileUpRejection();
  ~EpdPileUpRejection();
  
  void SetEpdPileUpCut(int index);  
  void SetTofPileUpCut(int index);  
  void PileUpEvent(StPicoDst * mPicoDst, Int_t & GoodEpd,Int_t & GoodTof);

  Float_t EpdNMip(StPicoDst * mPicoDst);
  Int_t Centrality(int val);
  std::pair<int,int> FxtMultAndPiPDuCount(StPicoDst * picoDst);

 private:
    
  StPicoDstMaker * mPicoDstMaker;
  StEpdGeom * mEpdGeom;
  StEpdEpInfo *mEpdEpInfo;

  std::vector< std::vector< float > > epdCutArray;
  std::vector< std::vector< float > > tofCutArray;

  StPicoDst * _mPicoDst;

  void FillArrays();

  float  getTofBeta   (StPicoDst * mPicoDst, StPicoTrack const* const trk, TVector3 const& vtxTV3) const;
  bool   isTofPion    (StPicoTrack const* const trk, float beta, TVector3 const& vtx) const;
  bool   isTofProton  (StPicoTrack const* const trk, float beta, TVector3 const& vtx) const;
  bool   isTofDeuteron(StPicoTrack const* const trk, float beta, TVector3 const& vtx) const;


  ClassDef(EpdPileUpRejection, 1);
  
};

#endif
