/*!
 * \class  StMuAnalysisMaker
 * \brief  A typical Analysis Class for MuDst
 * \author Wei-Ming Zhang, KSU, Mar 2004
 *
 * This is an example of a maker to perform analysis using MuDst.
 *
 * $Id: StMuAnalysisMaker.cxx,v 1.1 2004/08/10 16:09:11 perev Exp $
 * -------------------------------------------------------------------------
 * $Log: StMuAnalysisMaker.cxx,v $
 * Revision 1.1  2004/08/10 16:09:11  perev
 * new GridCollector stuff
 *
 * -------------------------------------------------------------------------
 */
//
//  Include header files.
#include "TFile.h"
#include "StMessMgr.h"
#include "TH1.h"
#include "TH2.h"
#include <map>
#include <assert.h>
#include <algorithm>

#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDebug.h"
#include "StBTofHeader.h"
#include "StMuAnalysisMaker.h"

ClassImp(StMuAnalysisMaker)

// The constructor. Initialize data members here.
StMuAnalysisMaker::StMuAnalysisMaker(const Char_t *name) : StMaker(name)
{
   mEventCounter = 0;
   mFile = 0;
}

StMuAnalysisMaker::~StMuAnalysisMaker()
{
   /* noop */
}

//  Called once at the beginning.
Int_t StMuAnalysisMaker::Init()
{
   //  Output file and histogram booking
  if(mFileName == NULL) mFileName = "muAnalysis.root";

   const char* EventsVarlist = "runid:eventid:" //
                               "vtx:vty:vtz:vpdVz:"
                               "gref:ref:nGlTracks:nPrTracks:nGlTracks2:nPrTracks2:nRecoTracks:nMCTracks:nDataTracks:dataCuts:dataDCA:FXTMult3:bField";
   EventsNtuple = new TNtuple("Events", "Events", EventsVarlist);

   //===================================================================
   // Nest is for Embedding Efficiency and momentum resolution
   //===================================================================
   const char* EmbeddingTracksVarlist = "runid:eventid:gref:ref:nPrimary:nReco:nData:nMC:FXTMult3:vtx:vty:vtz:vpdVz:ispion:iskaon:isproton:" //event vertex and flag of track pi/k
                                        "mcPt:mcP:mcEta:mcY:mcPhi:mcGeantId:mcHitsTpc:mcIdVx:" // MC track
                                        "rcId:rcPt:rcP:rcEta:rcPhi:rcCharge:rcNfit:rcNposs:rcNdedx:rcDedx:rcNsigPi:rcNsigK:rcDca:rcDcaXy:rcDcaZ:map0:map1:map2"; //  RC track
   EmbeddingTracksNtuple = new TNtuple("EmbeddingTracks", "EmbeddingTracks", EmbeddingTracksVarlist);

   return StMaker::Init();
}

//  Called every event after Make().
void StMuAnalysisMaker::Clear(Option_t *opt)
{
   StMaker::Clear();
}

//  Called once at the end.
Int_t StMuAnalysisMaker::Finish()
{
//  Summarize the run.
   cout << "StMuAnalysisMaker::Finish()\n";
   cout << "\tProcessed " << mEventCounter << " events." << endl;

//  Output histograms
   cout <<  "file name :" << mFileName.c_str() << endl;
   mFile =  new TFile(mFileName.c_str(), "RECREATE");
   cout << "\tHistograms will be stored in file '"
        <<  mFileName.c_str() << "'" << endl;

   EventsNtuple->Write();
   EmbeddingTracksNtuple->Write();

//  Write histos to file and close it.
   if (mFile)
   {
      mFile->Write();
      mFile->Close();
   }

   return kStOK;
}

//  This method is called every event.
Int_t StMuAnalysisMaker::Make()
{
   mEventCounter++;  // increase counter

//    DEBUGVALUE2(mEventCounter);
//  Get MuDst
   StMuDst* mu;
   mu = (StMuDst*) GetInputDS("MuDst");
//    DEBUGVALUE2(mu);

   if (!mu)
   {
      gMessMgr->Warning() << "StMuAnalysisMaker::Make : No MuDst" << endm;
      return kStOK;        // if no event, we're done
   }

//  Check StMuEvent branch
   StMuEvent* muEvent;
   muEvent = (StMuEvent*) mu->event();

   if(! muEvent->triggerIdCollection().nominal().isTrigger(620052)) return kStOK; 

   if (muEvent)
   {
      int refMult = muEvent->refMult();
      // cout << "refMult: " << refMult << endl;
      //      mRefMult->Fill(refMult);
   }

//  Printout information of StMuEvent
   muEventInfo(*mu);

//  MC track information, fix the mapping between MC and RC track // from MuTrack to find the MC Track
   TClonesArray *mcTracks = mu->mcArray(MCTrack);
   int nMcTracks = mcTracks->GetEntriesFast();

   map<int, int> index2McTrack;

   // printf("Mc track # = %d\n", nMcTracks);
   for (int i = 0; i < nMcTracks; i++)
   {
      StMuMcTrack *mcT = (StMuMcTrack *)mcTracks->UncheckedAt(i);
      if (!mcT) continue;
      int id = mcT->Id();
      int gId = mcT->GePid();
      //if(gId!=8) continue; // not a pion
      int idVtx = mcT->IdVx();
      //if(idVtx!=1) continue;  // not from primary vertex
      index2McTrack.insert(pair<int, int>(id, i));
   }

   int nRecoTracks = 0;
   int nDataTracks = 0;
   int dataCuts = 0;
   int dataDCA = 0;
   int FXTMult3 = 0;

   // build mapping between RC and MC track // from Mc Track to find the RC Track// the returned tracks can be more than 1
   multimap<int, int> index2RcTrack;


   // printf("global track # = %d\n", mu->globalTracks()->GetEntries());
   for (int i = 0; i < mu->globalTracks()->GetEntries(); i++)
   {
      StMuTrack *muTrack = (StMuTrack*) mu->globalTracks(i);
      if (!muTrack) continue;
      int idTruth = muTrack->idTruth();
      int qaTruth = muTrack->qaTruth();
      // if(qaTruth>0) cout << " i = " << std::setw(6) << i << " qaTruth = " << std::setw(6) << qaTruth << endl;
      if (idTruth <= 0) continue;
      if (idTruth > 10000) continue; // reconstructed tracks
      // cout << std::setw(6) << i << " " << std::setw(6) << muTrack->pt() << " :: " <<  std::setw(6) << idTruth-1 << " "<< std::setw(6)  << endl;
      index2RcTrack.insert(pair<int, int>(idTruth - 1, i));
   }


   //-----------------//
   // build mapping between RC-Prim and MC track // from Mc Track to find the RC Prim Track// the returned tracks can be more than 1
   multimap<int, int> index2PrimRcTrack;

   // printf("primary track # = %d\n", mu->primaryTracks()->GetEntries());
   for (int i = 0; i < mu->primaryTracks()->GetEntries(); i++)
   {
      StMuTrack *muTrack = (StMuTrack*) mu->primaryTracks(i);
      if (!muTrack) continue;
      int idTruth = muTrack->idTruth();
      int qaTruth = muTrack->qaTruth();
      // if(qaTruth>0) cout << " i = " << std::setw(6) << i << " qaTruth = " << std::setw(6) << qaTruth << endl;
      if (idTruth <= 0){
	nDataTracks++;
	if( muTrack->dca().mag() <= 3.0) dataDCA++;
	if( muTrack->dca().mag() <= 3.0 && (1.0*muTrack->nHitsFit() / muTrack->nHitsPoss()) >= 0.51) dataCuts++;
	if( muTrack->dca().mag() <= 3.0 && (1.0*muTrack->nHitsFit() / muTrack->nHitsPoss()) >= 0.51 && muTrack->nSigmaProton() <= -3.0) FXTMult3++;
	continue;
      }
      if (idTruth > 10000){
	nDataTracks++;
	if( muTrack->dca().mag() <= 3.0) dataDCA++;
	if( muTrack->dca().mag() <= 3.0 && (1.0*muTrack->nHitsFit() / muTrack->nHitsPoss()) >= 0.51) dataCuts++;
	if( muTrack->dca().mag() <= 3.0 && (1.0*muTrack->nHitsFit() / muTrack->nHitsPoss()) >= 0.51 && muTrack->nSigmaProton() <= -3.0) FXTMult3++;
	continue; // reconstructed tracks
      }
      // cout << std::setw(6) << i << " " << std::setw(6) << muTrack->pt() << " :: " <<  std::setw(6) << idTruth-1 << " "<< std::setw(6)  << endl;
      nRecoTracks++;
      index2PrimRcTrack.insert(pair<int, int>(idTruth - 1, i));
   }


   // =============  Build map between global and primary tracks from proper vertex
   TClonesArray *PrimaryVertices   = mu->primaryVertices();
   map<int, int> indexGl2Pr;

   for (int k = 0; k < mu->primaryTracks()->GetEntries(); k++)
   {
      StMuTrack *pTrack = (StMuTrack *) mu->primaryTracks(k);
      if (!pTrack) continue;
      int l = pTrack->vertexIndex();
      if (l < 0) continue;
      StMuPrimaryVertex *Vtx = (StMuPrimaryVertex *) PrimaryVertices->UncheckedAt(l);
      if (Vtx->idTruth() != 1) continue;
      int kg = pTrack->index2Global();

      indexGl2Pr.insert(pair<int, int>(kg, k));
   }



   //  Check global track branches
   StMuTrack* muTrack;
   int nTracks;
   nTracks = mu->globalTracks()->GetEntries();
   // printf("Global track # = %d\n", nTracks);
   for (int l = 0; l < nTracks; l++)
   {
      muTrack = (StMuTrack*) mu->globalTracks(l);
      if (!muTrack) continue;
      /*cout << "muGloTrack.flag : "  << muTrack->flag() << endl;  if (muTrack->idTruth())*/
      // { cout << "  idTruth : " << std::setw(6) << muTrack->idTruth() << endl; }
      if (!accept(muTrack)) continue;
      //      mGlobalPt->Fill(muTrack->pt());
      const StThreeVectorF mom = muTrack->momentum();
      const float pt = mom.perp();
      const float eta = mom.pseudoRapidity();
      // if (pt < 0.15 || fabs(eta) > 2.0) continue;
      if (fabs(eta) > 2.0) continue;

      int k = indexGl2Pr[l];
      float ppt = -1;
      if (k >= 0 && mu->primaryTracks()->GetEntries() >0)
      {
         StMuTrack *pTrack = (StMuTrack *) mu->primaryTracks(k);
         ppt = pTrack->pt();
      }

      unsigned long long map0 = muTrack->topologyMap().data(0);
      unsigned long long map1 = muTrack->topologyMap().data(1);
      unsigned long long map2 = muTrack->topologyMap().data(2);
   
      int idTruth = muTrack->idTruth();
      if (idTruth < 0) continue;
      if (idTruth > 10000) // reconstructed tracks
      {
	//         mRcTPCPtEta->Fill(pt, eta);
      }
      else
      {
         int index2Mc = index2McTrack[idTruth];
         if (index2Mc >= 0)
         {
            StMuMcTrack *mcT = (StMuMcTrack *)mcTracks->UncheckedAt(index2Mc);
            if (mcT)
            {
               float pt_mc = mcT->pT();
	       //	       mPtCorr->Fill(pt_mc, pt - pt_mc);
	       //               mMcTPCPtEta->Fill(pt, eta);
            }
         } // end if (index2Mc)
      }
   }



   //  Check primary track branches
   StMuTrack* pTrack;
   int nPTracks;
   nPTracks = mu->primaryTracks()->GetEntries();
   // printf("Primary track # = %d\n", nPTracks);
   for (int l = 0; l < nPTracks; l++)
     {
       pTrack = (StMuTrack*) mu->primaryTracks(l);
       if (!pTrack) continue;
       // cout << "muPriTrack.flag : "  << pTrack->flag() << endl;
       //int kg = pTrack->index2Global(); // find the correspond global track
       if (!accept(pTrack)) continue;
       //       mPrimaryPt->Fill(pTrack->pt());
       const StThreeVectorF mom = pTrack->momentum();
       const float pt = mom.perp();
       const float eta = mom.pseudoRapidity();
       // if (pt < 0.15 || fabs(eta) > 2.0) continue;
       if (fabs(eta) > 2.0) continue;
       
       int idTruth = pTrack->idTruth();
       if (idTruth < 0) continue;
       if (idTruth > 10000) continue; // reconstructed tracks
       
       int index2Mc = index2McTrack[idTruth];
       if (index2Mc >= 0)
	 {
	   StMuMcTrack *mcT = (StMuMcTrack *)mcTracks->UncheckedAt(index2Mc);
	   if (mcT)
	     {
	       float pt_mc = mcT->pT();
	     }
	 }
       
     }
   

   
   //  Check MC track branches
   for (int i = 0; i < nMcTracks; i++)
     {
       StMuMcTrack *mcT = (StMuMcTrack*)mcTracks->UncheckedAt(i);
       if (!mcT) continue;
       int id = mcT->Id() - 1;
       pair<multimap<int, int>::iterator, multimap<int, int>::iterator> ret;
       // ret = index2RcTrack.equal_range(id);
       ret = index2PrimRcTrack.equal_range(id);
       multimap<int, int>::iterator it;
       int kg = -1;
       int count = 0;
       float deltaPt = 999;//intial a large num
       int Tmpkg = -1;
       int TmpqaTruth = -1;
       StMuTrack *muTrack;
       // for (it = index2RcTrack.equal_range(id).first; it != index2RcTrack.equal_range(id).second; ++it, ++count)
       for (it = index2PrimRcTrack.equal_range(id).first; it != index2PrimRcTrack.equal_range(id).second; ++it, ++count)
      {
         kg = (*it).second;
         // muTrack = (StMuTrack*) mu->globalTracks(kg);
         muTrack = (StMuTrack*) mu->primaryTracks(kg);
         if (!muTrack) continue;
         // if (fabs(muTrack->pt() - mcT->pT()) < deltaPt) //use delta pt select rc track
         if (muTrack->qaTruth() > TmpqaTruth) //use delta pt select rc track
         {
            deltaPt = fabs(muTrack->pt() - mcT->pT());
            TmpqaTruth = muTrack->qaTruth();
            Tmpkg = kg;
         }
      }
       if (count > 0) //clone tracks  //   // count > 0 means one mc track have multiple global tracks, here use the closest pT rc track as associated one
	 {
	   // muTrack = (StMuTrack*) mu->globalTracks(Tmpkg);
	   muTrack = (StMuTrack*) mu->primaryTracks(Tmpkg);
	 }
       
       if (Tmpkg > -1) //Tmpkg must be large than -1. then it means the Rc track was found
	 {

      }
   }



   //========================================================================
   //
   //Next is for some basic Event information NTuple
   //
   //========================================================================
   // fill ntuple
   float array2[50];
   int idx2 = 0;
   array2[idx2++] = muEvent->runId();
   array2[idx2++] = muEvent->eventId();
   // Event
   array2[idx2++] = muEvent->primaryVertexPosition().x();//vx
   array2[idx2++] = muEvent->primaryVertexPosition().y();//vy
   array2[idx2++] = muEvent->primaryVertexPosition().z();//vz
   array2[idx2++] = muEvent->vpdVz();//vpd vz
   //Refit Vertes //
   array2[idx2++] = muEvent->grefmult();
   array2[idx2++] = muEvent->refMult();
   array2[idx2++] = mu->globalTracks()->GetEntries();
   array2[idx2++] = mu->primaryTracks()->GetEntries();
   array2[idx2++] = mu->array(muGlobal)->GetEntriesFast();
   array2[idx2++] = mu->array(muPrimary)->GetEntriesFast();
   array2[idx2++] = nRecoTracks;
   array2[idx2++] = nMcTracks;
   array2[idx2++] = nDataTracks;
   array2[idx2++] = dataCuts;
   array2[idx2++] = dataDCA;
   array2[idx2++] = FXTMult3;
   array2[idx2++] = muEvent->magneticField();

   EventsNtuple->Fill(array2);

   if (!isGoodEvent(*mu)) return kStOk;


   //========================================================================
   //
   //Next is for Extract Embedding Efficiency for Singal Track (pi/k) and momentum resolution
   //
   //========================================================================
   for (int i_Mc = 0; i_Mc < nMcTracks; ++i_Mc)
   {
      StMuMcTrack *mcTrk = (StMuMcTrack *)mcTracks->UncheckedAt(i_Mc);
      if (!mcTrk)  continue;

      if (!isGoodEmbeddingTrack(mcTrk)) continue;
      int idVtx = mcTrk->IdVx(); // cout << std::setw(6) << idVtx << endl;
      if (idVtx != 1) continue; // not from primary vertex

      bool isPion = isPionTrack(mcTrk) ;
      bool isKaon = isKaonTrack(mcTrk) ;
      bool isProton = isProtonTrack(mcTrk) ;

      if (!isPion && !isKaon && !isProton) continue;

      int id = mcTrk->Id() - 1;

      pair<multimap<int, int>::iterator, multimap<int, int>::iterator> ret;
      // ret = index2RcTrack.equal_range(id);
      ret = index2PrimRcTrack.equal_range(id);
      multimap<int, int>::iterator it;
      int kg = -1;
      int count = 0;
      float deltaPt = 999;//intial a large num
      int Tmpkg = -1;
      int TmpqaTruth = -1;
      StMuTrack *muTrack = 0;
      // for (it = index2RcTrack.equal_range(id).first; it != index2RcTrack.equal_range(id).second; ++it, ++count)
      for (it = index2PrimRcTrack.equal_range(id).first; it != index2PrimRcTrack.equal_range(id).second; ++it, ++count)
      {
         kg = (*it).second;
         // muTrack = (StMuTrack*) mu->globalTracks(kg);
         muTrack = (StMuTrack*) mu->primaryTracks(kg);
         if (!muTrack) continue;
         // if (fabs(muTrack->pt() - mcTrk->pT()) < deltaPt)
         if (muTrack->qaTruth() > TmpqaTruth) //use delta pt select rc track
         {
            deltaPt = fabs(muTrack->pt() - mcTrk->pT());
            TmpqaTruth = muTrack->qaTruth();
            Tmpkg = kg;
         }
      }
      if (count > 0) //clone tracks  //   // count > 0 means one mc track have multiple global tracks, here use the closest pT rc track as associated one
      {
         // muTrack = (StMuTrack*) mu->globalTracks(Tmpkg);
         muTrack = (StMuTrack*) mu->primaryTracks(Tmpkg);
      }

      if (Tmpkg > -1) //Tmpkg must be large than -1. then it means the Rc track was found// otherwise the asscoiate RC was not found
      {
         // cout << " i = " << std::setw(6) << i_Mc << " k = " << std::setw(6) << Tmpkg << " trkID = " << std::setw(6) << muTrack->id()
         //   << " pt = " << std::setw(6) << muTrack->pt() << " pt_mc = " << std::setw(6) << mcTrk->pT()  <<  endl; // count > 0 then one mc track have multiple global tracks
      }

      bool rcflag = (Tmpkg > -1); // cout << " i_Mc = " << std::setw(6) << i_Mc << " k = " << std::setw(6) << Tmpkg << " rcflag : " << std::setw(6) << rcflag << endl;


      // fill ntuple
      float array[200];
      int idx = 0;
      array[idx++] = muEvent->runId();
      array[idx++] = muEvent->eventId();
      array[idx++] = muEvent->grefmult();
      array[idx++] = muEvent->refMult();
      array[idx++] = mu->numberOfPrimaryTracks();
      array[idx++] = nRecoTracks;
      array[idx++] = nDataTracks;
      array[idx++] = nMcTracks;
      array[idx++] = FXTMult3;
      array[idx++] = muEvent->primaryVertexPosition().x();//vx
      array[idx++] = muEvent->primaryVertexPosition().y();//vy
      array[idx++] = muEvent->primaryVertexPosition().z();//vz
      array[idx++] = muEvent->vpdVz();//vpd vz
      // MC pion/kaon
      array[idx++] = isPion;
      array[idx++] = isKaon;
      array[idx++] = isProton;
      array[idx++] = mcTrk->pT();//pt
      array[idx++] = mcTrk->Ptot();//p
      array[idx++] = mcTrk->Eta();
      array[idx++] = mcTrk->Rapidity();
      array[idx++] = atan2(mcTrk->Pxyz().y(), mcTrk->Pxyz().x());
      array[idx++] = mcTrk->GePid();
      array[idx++] = mcTrk->NoHits();
      array[idx++] = mcTrk->IdVx();
      // RC pion/kaon
      array[idx++] = rcflag ? muTrack->id() : -9999.; //rc track Id for check
      // cout << " trkID = " << std::setw(6) << (rcflag ? muTrack->id() : -9999.)  << " muTrack? = " << std::setw(6) << muTrack << endl;
      fill_rcTrack(array, idx, mu, muEvent, muTrack, mcTrk);

      EmbeddingTracksNtuple->Fill(array);


   } // .. end mc tracks loop


   return kStOK;
}

//  A simple track filter
bool StMuAnalysisMaker::accept(StMuTrack* track)
{
//  check for positive flags.
   return track && track->flag() >= 0 && track->dcaGlobal().mag() < 3.0 && (int)track->nHitsFit() > 10; // && fabs(track->nSigmaPion())<2.0;
}

//  Prototype
void StMuAnalysisMaker::muEventInfo(const StMuDst& mu)
{
   StMuEvent* ev = mu.event();
   if (!ev) return;

   StThreeVectorF pVtx = ev->primaryVertexPosition();
   //   mVxVy->Fill(pVtx.x(), pVtx.y());
   float vzVpd = -999.;
   if (StBTofHeader* header = mu.btofHeader())
   {
      vzVpd = header->vpdVz();
   }
   // if (fabs(vzVpd) < 210.)
   {
     //      mVzCorr->Fill(pVtx.z(), vzVpd);
   }

   int gRefMult = ev->grefmult();
   int refMult = ev->refMultNeg() + ev->refMultPos();
   //   mRefCorr->Fill(refMult, gRefMult);

}

//-----------------------------------------------------------------------------
bool StMuAnalysisMaker::isGoodEvent(const StMuDst& mu)
{
   float vx = mu.event()->primaryVertexPosition().x();//vz
   float vy = mu.event()->primaryVertexPosition().y();//vz
   float vz = mu.event()->primaryVertexPosition().z();//vz
   float vpdVz = mu.event()->vpdVz();//vpdVz
   return ( vz<202 && vz>198 && sqrt(pow(vx,2.0)+pow(vy+2.0,2.0))<1.5);
   //   return (fabs(vz-200.) <= 2.0 && sqrt(pow(vx, 2.) + pow(vy + 2, 2.)) <= 2.0 );
}
//-----------------------------------------------------------------------------
bool StMuAnalysisMaker::isGoodEmbeddingTrack(StMuMcTrack* mcTrk)
{
   // return (mcTrk->pT() > 0.15 && std::fabs(mcTrk->Eta()) <= 2.0);
   return (std::fabs(mcTrk->Eta()) <= 2.0);
}
//-----------------------------------------------------------------------------
bool StMuAnalysisMaker::isPhiTrack(StMuMcTrack* mcTrk)
{
   if (!mcTrk) return false;
   int geantid = mcTrk->GePid();
   return (geantid == 34);   // phi + 151?
}
//-----------------------------------------------------------------------------
bool StMuAnalysisMaker::isPionTrack(StMuMcTrack* mcTrk)
{
   if (!mcTrk) return false;
   int geantid = mcTrk->GePid();
   return (geantid == 8  || geantid == 9);  // pion + -
   // return (geantid == 8 );  // pion +
}
//-----------------------------------------------------------------------------
bool StMuAnalysisMaker::isKaonTrack(StMuMcTrack* mcTrk)
{
   if (!mcTrk) return false;
   int geantid = mcTrk->GePid();
   return (geantid == 11  || geantid == 12) ;  // kaon + -
   // return (geantid == 12) ;  // kaon -
}
//---------------------------
bool StMuAnalysisMaker::isProtonTrack(StMuMcTrack* mcTrk)
{
   if (!mcTrk) return false;
   int geantid = mcTrk->GePid();
   return (geantid == 14  || geantid == 15);  // proton + -
   // return (geantid == 14 );  // proton +
}
//-----------------------------------------------------------------------------

void StMuAnalysisMaker::fill_rcTrack(float* array, int& idx, StMuDst* mu, StMuEvent* muEvent, StMuTrack* muTrack, StMuMcTrack* mcTrk)
{

   StThreeVectorF pVtx(-999., -999., -999.);
   pVtx = muEvent->primaryVertexPosition();
   StThreeVectorF  momentum(-999., -999., -999);
   if (muTrack) momentum = muTrack->momentum();

   array[idx++] = muTrack ? momentum.perp() : -9999.;
   array[idx++] = muTrack ? momentum.mag() : -9999.;
   array[idx++] = muTrack ? momentum.pseudoRapidity() : -9999.;
   array[idx++] = muTrack ? momentum.phi() : -9999.;
   array[idx++] = muTrack ? muTrack->charge() : -9999.;
   array[idx++] = muTrack ? muTrack->nHitsFit() : -9999.;
   array[idx++] = muTrack ? muTrack->nHitsPoss() : -9999.;
   array[idx++] = muTrack ? muTrack->nHitsDedx() : -9999.;

   // dedx info
   array[idx++] = muTrack ? muTrack->dEdx() : -9999.;
   array[idx++] = muTrack ? muTrack->nSigmaPion() : -9999.;
   array[idx++] = muTrack ? muTrack->nSigmaKaon() : -9999.;

   // dca
   // array[idx++] = muTrack ? muTrack->helix().geometricSignedDistance(pVtx) : -9999.;//dca
   // array[idx++] = muTrack ? muTrack->helix().geometricSignedDistance(pVtx.x(), pVtx.y()) : -9999.;//dcaXy
   // array[idx++] = muTrack ? (muTrack->helix().at(muTrack->helix().pathLength(pVtx.x(), pVtx.y()))).z() - pVtx.z() : -9999.;//dcaZ
   StMuTrack* gTrack = 0;
   if(muTrack) gTrack = (StMuTrack*) mu->globalTracks(muTrack->index2Global());
   array[idx++] = gTrack ? gTrack->helix().geometricSignedDistance(pVtx) : -9999.;//dca
   array[idx++] = gTrack ? gTrack->helix().geometricSignedDistance(pVtx.x(), pVtx.y()) : -9999.;//dcaXy
   array[idx++] = gTrack ? (gTrack->helix().at(gTrack->helix().pathLength(pVtx.x(), pVtx.y()))).z() - pVtx.z() : -9999.;//dcaZ

   //topoMap
   // unsigned long long map0 = muTrack->topologyMap().data(0);
   // unsigned long long map1 = muTrack->topologyMap().data(1);
   // unsigned long long map2 = muTrack->topologyMap().data(2);
   array[idx++] = muTrack ? (unsigned long long)muTrack->topologyMap().data(0) : -9999.; //dca
   array[idx++] = muTrack ? (unsigned long long)muTrack->topologyMap().data(1) : -9999.;//dca
   array[idx++] = muTrack ? (unsigned long long)muTrack->topologyMap().data(2) : -9999.;//dca

}


void StMuAnalysisMaker::setPart(string partCStr){

  TString partStr = Form("%s",partCStr.c_str());

  string pidSym;

  if( partStr.EqualTo("pip",TString::kIgnoreCase) ){
    part = 0;
    ch = 1;
    pidSym = "#pi^{+}";
    mass = 0.13957;
  }else if( partStr.EqualTo("pim",TString::kIgnoreCase) ){
    part = 0;
    ch = -1;
    pidSym = "#pi^{-}";
    mass = 0.13957;
  }else if( partStr.EqualTo("Kp",TString::kIgnoreCase) ){
    part = 1;
    ch = 1;
    pidSym = "K^{+}";
    mass = 0.49367;
  }else if( partStr.EqualTo("Km",TString::kIgnoreCase) ){
    part = 1;
    ch = -1;
    pidSym = "K^{-}";
    mass = 0.49367;
  }else if( partStr.EqualTo("pro",TString::kIgnoreCase) ){
    part = 2;
    ch = 1;
    pidSym = "p";
    mass = 0.93827;
  }

  pidSymbol = pidSym.c_str();

  return;

}
