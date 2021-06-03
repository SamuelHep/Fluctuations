/*!
 * \class  StMuAnalysisMaker
 * \brief  A typical Analysis Class for MuDst
 * \author Wei-Ming Zhang, KSU, Mar 2004
 *
 * This is an example of a maker to perform analysis using MuDst.
 *
 * $Id: StMuAnalysisMaker.h,v 1.2 2014/08/06 11:43:31 jeromel Exp $
 *
 * -------------------------------------------------------------------------
 * $Log: StMuAnalysisMaker.h,v $
 * Revision 1.2  2014/08/06 11:43:31  jeromel
 * Suffix on literals need to be space (later gcc compiler makes it an error) - first wave of fixes
 *
 * Revision 1.1  2004/08/10 16:09:11  perev
 * new GridCollector stuff
 *
 * -------------------------------------------------------------------------
 */
#ifndef StMuAnalysisMaker_hh
#define StMuAnalysisMaker_hh
//
//  Include files
#include "StMaker.h"
#include <string>
#include <vector>
#include "TNtuple.h"

//
//  Forward declarations
class StMuTrack;
class TFile;
class TH1D;
class TH2F;
class TH3F;
class StMuDst;
class StMuEvent;
class StMuTrack;
class StMuMcTrack;
class TNtuple;

#ifndef ST_NO_NAMESPACES
using std::string;
#endif
//
//  The class declaration. It innherits from StMaker.
class StMuAnalysisMaker : public StMaker
{
public:

  const int nCentBins = 9;
  //string pidSym = "#pi^{+}";
  
  int part = -999;
  int ch = -999;
  
  double mass;

  const char* pidSymbol;// = pidSym.c_str();

   StMuAnalysisMaker(const Char_t *name = "muAnalysis"); // constructor
   ~StMuAnalysisMaker();                                 // destructor

   void Clear(Option_t *option = ""); // called after every event to cleanup
   Int_t  Init();                   // called once at the beginning of your job
   Int_t  Make();                   // invoked for every event
   Int_t  Finish();                 // called once at the end
   void muEventInfo(const StMuDst &mu);

   void setOutputName(string a_name){mFileName = a_name; };

   void setPart(string partCStr);

   bool isGoodEvent(const StMuDst& mu);
   bool isGoodEmbeddingTrack(StMuMcTrack* mcTrk);
   bool isPhiTrack(StMuMcTrack* mcTrk);
   bool isPionTrack(StMuMcTrack* mcTrk);
   bool isKaonTrack(StMuMcTrack* mcTrk);
   bool isProtonTrack(StMuMcTrack* mcTrk);
   void fill_rcTrack(float* array, int& idx, StMuDst* mu, StMuEvent* muEvent, StMuTrack* muTrack, StMuMcTrack* mcTrk);


   virtual const char *GetCVS() const
   {
      static const char cvs[] = "Tag $Name:  $ $Id: StMuAnalysisMaker.h,v 1.2 2014/08/06 11:43:31 jeromel Exp $ built " __DATE__ " " __TIME__ ;
      return cvs;
   }

private:

// data member
   int        mEventCounter;  //!
   string     mFileName;      //!
   TFile      *mFile;         //!

   int nRapidityBins;
   double rapidityMin;
   double rapidityMax;

   TNtuple* EventsNtuple;
   TNtuple* EmbeddingTracksNtuple;

// method (a simple track filter)
   bool accept(StMuTrack*);            // and this is used to select tracks

   ClassDef(StMuAnalysisMaker, 0)
};
#endif
