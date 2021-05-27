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
class TH1D;
class TF1;

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

  int LoadTPCEff(TString filename);
  int GetPionEffBin( double y );
  double GetPionEfficiency(double y, double pt);
  Double_t ParticleRapidity( Double_t pt, Double_t pz, Double_t mass );

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

  std::vector<TH1D*> _h1d_eff;
  std::vector<TF1*> _tf1_eff;


  ClassDef(StFemtoDstMaker,1)

};


#endif
