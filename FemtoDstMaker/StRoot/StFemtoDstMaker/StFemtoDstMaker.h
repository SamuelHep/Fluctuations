#ifndef ST_FEMTO_DSTMAKER_h
#define ST_FEMTO_DSTMAKER_h

// C++ headers
#include <vector>

// ROOT headers
#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "StPicoEvent/StPicoTrack.h"

#include "StRoot/StFemtoTrack/StFemtoTrack.h"
#include "StRoot/StFemtoEvent/StFemtoEvent.h"

#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/EpdPileUpRejection/EpdPileUpRejection.h"
#include "TVector3.h"

class StFemtoTrack;

//_________________
class StFemtoDstMaker : public StMaker {

 public:
  /// Default constructor
  StFemtoDstMaker(char * name);

  /// Destructor
  virtual ~StFemtoDstMaker();
  /// Print some event information

  Int_t Init();
  Int_t Finish();
  Int_t Make();

  void InitRunIndex();
  int GetRunIndex(int runId);
  bool GoodRun(int runId); 

  void SetRunFile(TString val){runNumFile = val;};
  void SetFileIndex(char *val) {mFileIndex=val;}
  void SetOutDir(char *val) {mOutDir=val;}

  Float_t GetT0(StPicoTrack const* const trk, TVector3 const& vtxTV3,float mass) const;

 private:
  
  StPicoDstMaker * mPicoDstMaker;
  EpdPileUpRejection * mEpdPileUpRejection;

  TString runNumFile;
  std::vector <int> run_vec;

  char *mFileIndex;
  char *mOutDir;
  StFemtoEvent * mFemtoEvent;
  TTree * fDstTree;
  TFile * mOutfile;

  int fullEventCounter;

  ClassDef(StFemtoDstMaker,1)

};


#endif
