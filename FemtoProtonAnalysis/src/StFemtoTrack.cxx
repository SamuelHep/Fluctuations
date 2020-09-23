
#include "TObject.h"
#include <iostream>
#include "StFemtoTrack.h"

ClassImp(StFemtoTrack)

//___________________________________________________
StFemtoTrack::StFemtoTrack() : TObject(),
  mPx(-999), mPy(-999), mPz(-999), mNHitsFit(-122), mDcaX(-999),mDcaY(-999),mDcaZ(-999), 
			       mNSigmaProton(-122), mNSigmaPion(-122), mTofMass(-999), mTofT0(-999)
{
}

//___________________________________________________
StFemtoTrack::StFemtoTrack(const StFemtoTrack &track) : TObject()
{

  mPx = track.mPx;
  mPy = track.mPy;
  mPz = track.mPz;
  
  mNHitsFit = track.mNHitsFit;

  mDcaX = track.mDcaX;
  mDcaY = track.mDcaY;
  mDcaZ  = track.mDcaZ;
  
  mNSigmaProton = track.mNSigmaProton;
  mNSigmaPion = track.mNSigmaPion;
  mTofMass = track.mTofMass;
  mTofT0 = track.mTofT0;

}


void StFemtoTrack::PrintTrack()
{

  std::cout << "~TRACK~" << std::endl;
  std::cout << "   Px = " << mPx/1000.0 << std::endl;
  std::cout << "   Py = " <<  mPy/1000.0 << std::endl;
  std::cout << "   Pz = " <<  mPz/1000.0 << std::endl;
  
  std::cout << "   NHitsFit = " <<  ((int) mNHitsFit) << std::endl;

  std::cout << "   DcaX  = " <<  mDcaX/1000.0  << std::endl;
  std::cout << "   DcaY = " <<  mDcaY/1000.0 << std::endl;
  std::cout << "   DcaZ = " <<  mDcaZ/1000.0 << std::endl;
  
  std::cout << "   NSigmaProton = " <<  mNSigmaProton/100.0 << std::endl;
  std::cout << "   NSigmaPion = " <<  mNSigmaPion/100.0 << std::endl;
  std::cout << "   TofMass = " <<  GetTofMass() << std::endl;
  std::cout << "   TofT0 = " <<  GetTofT0() << std::endl;

}


