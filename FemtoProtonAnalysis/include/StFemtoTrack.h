/**
 * \class StFemtoTrack
 * \brief Holds information about track parameters
 *
 * The class stores information about the tracks reconstructed in TPC
 */

#ifndef StFemtoTrack_h
#define StFemtoTrack_h

// C++ headers
#include <cmath>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

// StFemtoTrack

//_________________
class StFemtoTrack : public TObject {

 public:
  /// Default constructor
  StFemtoTrack();

  /// Copy constructor
  StFemtoTrack(const StFemtoTrack &track);
  /// Destructor
  virtual ~StFemtoTrack(){};
  /// Print track parameters
  //  virtual void Print(const Char_t *option = "") const;
  void PrintTrack();


  //Getters
  Float_t GetPx(){return (Float_t) mPx/1000.0;};
  Float_t GetPy(){return (Float_t) mPy/1000.0;};
  Float_t GetPt(){return sqrt( pow(mPx/1000.0,2) + pow(mPy/1000.0,2) );};
  Float_t GetPz(){return (Float_t) mPz/1000.0;};
  Short_t GetNHitsFit(){return abs( (Short_t) mNHitsFit ) ;};
  Short_t GetCharge(){return  ((Short_t) mNHitsFit > 0) ? 1 : -1 ;};

  Float_t GetDcaX(){return (Float_t) mDcaX/1000.0;};
  Float_t GetDcaY(){return (Float_t) mDcaY/1000.0;};
  Float_t GetDcaZ(){return (Float_t) mDcaZ/1000.0;};


  Float_t GetDca(){return sqrt( pow( mDcaX/1000.0,2) + pow(mDcaY/1000.0,2) + pow( mDcaZ/1000.0,2) );};

  Float_t GetNSigmaProton(){return ((Float_t)mNSigmaProton)/100.0;};
  Float_t GetNSigmaPion(){return ((Float_t)mNSigmaPion)/100.0;};
  Float_t GetTofMass(){return (Float_t) (mTofMass>=0) ? mTofMass/1000.0 : -999;};
  Float_t GetTofT0(){return (Float_t) (mTofMass>=0) ? mTofT0/1000.0 : -999;};

  //Setters
  void SetPx(Float_t val){mPx = val*1000;};
  void SetPy(Float_t val){mPy = val*1000;};
  void SetPz(Float_t val){ mPz = val*1000;};
  void SetNHitsFitAndCharge(Float_t val1,Float_t val2){ mNHitsFit = val1*val2;};
  
  void SetDcaX(Float_t val){ mDcaX = val*1000;};
  void SetDcaY(Float_t val){ mDcaY = val*1000;};
  void SetDcaZ(Float_t val){ mDcaZ = val*1000;};

  void SetNSigmaProton(Float_t val){ mNSigmaProton= val*100; };
  void SetNSigmaPion(Float_t val){ mNSigmaPion = val*100; };
  void SetTofMass(Float_t val){ mTofMass = val*1000;};
  void SetTofT0(Float_t val){ mTofT0 = val*1000;};

 protected:

  Short_t mPx;
  Short_t mPy;
  Short_t mPz;
  Char_t  mNHitsFit;
  Short_t mDcaX;
  Short_t mDcaY;
  Short_t mDcaZ;
  Short_t  mNSigmaProton;
  Short_t  mNSigmaPion;
  Short_t mTofMass;
  Short_t mTofT0;


  ClassDef(StFemtoTrack,1)

};



#endif
