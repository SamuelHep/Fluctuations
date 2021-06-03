#include "TStopwatch.h"
#include<string>

class StMaker;
class StMuDstMaker;

class StChain;
StChain *chain=0;

void loadLibs();


void doMuDst(const Char_t *inputFileList = "test.list", const Char_t *outputFileName = "testing.root", Char_t *partStr="pro", int nEvents = 10000000)
{
   
  TStopwatch*   stopWatch = new TStopwatch();
  stopWatch->Start();

  gROOT->Macro("loadMuDst.C");
  //  loadLibs();
  gSystem->Load("StMuAnalysisMaker");
  //gSystem->Load("StMuDstMaker");

  chain = new StChain("StChain");
  chain->SetDebug();

  StMuDstMaker* MuDstMaker = new StMuDstMaker(0, 0, "", inputFileList, "MuDst", 10000);

  MuDstMaker->SetStatus("*", 1);
   
  StMuAnalysisMaker *anaMaker = new StMuAnalysisMaker();

  anaMaker->setOutputName(outputFileName);
  anaMaker->setPart(partStr);

  Int_t initStat = chain->Init();
  if (initStat) chain->Fatal(initStat, "Failure During Init()");

  cout << "chain->Init();" << endl;

  int total = MuDstMaker->chain()->GetEntries();
  //  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i = 0; i < nEvents; i++)
  {
    if (i % 100 == 0)
      cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret)
    {
      cout << "Bad return code!" << iret << endl;
      break;
    }

    //    total++;
  }

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!" << endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  delete chain;
  
  stopWatch->Stop();
  stopWatch->Print();

}



void loadLibs()
{
  //  gSystem->Load("libTable");
  //  gSystem->Load("libPhysics");
  gSystem->Load("St_base"); //
  gSystem->Load("StChain"); //
  gSystem->Load("St_Tables"); //
  gSystem->Load("StUtilities"); //        // new addition 22jul99
    gSystem->Load("StTreeMaker"); //
    gSystem->Load("StIOMaker"); //
    gSystem->Load("StarClassLibrary"); //
  gSystem->Load("StTriggerDataMaker"); // new starting from April 2003
  gSystem->Load("StBichsel"); //
  gSystem->Load("StEvent"); //
  gSystem->Load("StEventUtilities");  //
  gSystem->Load("StDbLib"); //
  gSystem->Load("StEmcUtil"); //
  gSystem->Load("StTofUtil"); //
  //  gSystem->Load("StPmdUtil");
  //  gSystem->Load("StPreEclMaker");
  gSystem->Load("StStrangeMuDstMaker"); //
  gSystem->Load("StMuDSTMaker"); //
  gSystem->Load("libStarAgmlUtil");

  /*
  gSystem->Load("StTpcDb");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("libgen_Tables");
  gSystem->Load("libsim_Tables");
  gSystem->Load("libglobal_Tables");
  gSystem->Load("StEmcTriggerMaker");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StTriggerUtilities");

  gSystem->Load("StMagF");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");
  */
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  gSystem->ListLibraries();
}
