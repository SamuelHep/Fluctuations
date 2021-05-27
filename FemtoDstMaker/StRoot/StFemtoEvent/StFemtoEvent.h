#ifndef StFemtoEvent_h
#define StFemtoEvent_h

// C++ headers
#include <vector>
#include <iostream>
// ROOT headers
#include "TObject.h"
#include "TClonesArray.h"
#include "StRoot/StFemtoTrack/StFemtoTrack.h"
#include "TVector3.h"


class StFemtoTrack;
//_________________
class StFemtoEvent : public TObject {

 public:
  /// Default constructor
  StFemtoEvent();
  /// Copy constructor
  StFemtoEvent(const StFemtoEvent &event);
  /// Destructor
  virtual ~StFemtoEvent();
  /// Print some event information

  //  virtual void Print(const Char_t *option = "") const;
  void PrintEvent();
  void ClearEvent();

  //
  // Getters
  //

  Int_t GetRunID(){return mRunID;};
  Short_t GetFxtMult(){return mFxtMult;};
  Float_t GetFxtMultTofMatch(){return mFxtMultTofMatch;};
  Short_t GetFxtMult3(){return mFxtMult3;};

  Float_t GetVz(){return (Float_t) mVz;};
  Float_t GetVxy(){return (Float_t) mVxy/1000.0;};
  
  Float_t GetNMip(){ return mNMip; };
  Short_t GetPiPDu(){ return mPiPDu; };

  //  Bool_t IsFullEvent(){ return mFullEvent;};

  //
  // Setters
  //

  void SetRunID(Int_t val){mRunID = val;};
  void SetFxtMult(Short_t val){mFxtMult = val;};
  void SetFxtMultTofMatch(Float_t val){mFxtMultTofMatch =val;};
  void SetFxtMult3(Short_t val){mFxtMult3 = val;};

  void SetVz(Float_t val){ mVz = val;};
  void SetVxy(Float_t val){ mVxy = val*1000;};
  void SetVr(Float_t val){ mVxy = val*1000;};
  
  void SetEpdNMip(Float_t val){ mNMip = val;};
  void SetPiPDu(int val){ mPiPDu = val;};
  //  void SetFullEvent(Bool_t val){ mFullEvent = val;};
  
  Int_t GetEntries(){return mFemtoTrackArray.size();}
  StFemtoTrack GetFemtoTrack(int i){return (StFemtoTrack) mFemtoTrackArray[i];}
  void SetStFemtoTrackArray(std::vector< StFemtoTrack >  val){mFemtoTrackArray = val;};


 private:

  Int_t mRunID;

  Short_t mFxtMult;
  Float_t mFxtMultTofMatch;
  Short_t mFxtMult3;

  Float_t mVz;
  Short_t mVxy;

  Float_t mNMip;
  Short_t mPiPDu;

  //  Bool_t mFullEvent;
  
  std::vector< StFemtoTrack>  mFemtoTrackArray;


  ClassDef(StFemtoEvent,1)

};

#endif
