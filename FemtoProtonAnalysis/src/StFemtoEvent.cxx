#include "StFemtoTrack.h"
#include "StFemtoEvent.h"


ClassImp(StFemtoEvent)

//____________________________________________________________
StFemtoEvent::StFemtoEvent(): TObject(),
  mRunID(0), mFxtMult(-999), mFxtMultTofMatch(-999), mFxtMult3(-999), mVz(-999), mVxy(-999),
  mNMip(-999), mPiPDu(-999) // mFullEvent(0) //mFemtoTrackArray(NULL)
{
}

//____________________________________________________________
StFemtoEvent::StFemtoEvent(const StFemtoEvent &event): TObject()
{
  mRunID = event.mRunID;
  mFxtMult = event.mFxtMult;
  mFxtMultTofMatch = event.mFxtMultTofMatch;
  mFxtMult3 = event.mFxtMult3;
  mVz = event.mVz;
  mVxy = event.mVxy;

  mNMip = event.mNMip;
  mPiPDu = event.mPiPDu;

  //  mFullEvent = event.mFullEvent;

  mFemtoTrackArray = event.mFemtoTrackArray;
}

//____________________________________________________________
StFemtoEvent::~StFemtoEvent(){}

void StFemtoEvent::ClearEvent()
{

  mRunID = 0;
  mFxtMult = -999;
  mFxtMultTofMatch = -999;
  mFxtMult3 = -999;
  mVz = -999;
  mVxy = -999;
  mNMip = -999;
  mPiPDu = -999;
  //  mFullEvent = 0;

  mFemtoTrackArray.clear();

}

void StFemtoEvent::PrintEvent()
{

  std::cout << "~EVENT~" << std::endl;
  std::cout << "   RunID    = " <<  mRunID << std::endl;
  std::cout << "   FxtMult  = " <<  mFxtMult << std::endl;
  std::cout << "   FxtMultTofMatch  = " <<  mFxtMultTofMatch << std::endl;
  std::cout << "   FxtMult3 = " <<  mFxtMult3 << std::endl;
  std::cout << "   Vz       = " <<  mVz << std::endl;
  std::cout << "   Vxy      = " <<  mVxy/1000.0 << std::endl;
  std::cout << "   NMip     = " <<  mNMip << std::endl;
  std::cout << "   PiPDu    = " <<  mPiPDu << std::endl;
  //  std::cout << "   FullEvent= " <<  mFullEvent << std::endl;
  //  std::cout << "   NTracks  = " <<  GetEntries() << std::endl;

}
