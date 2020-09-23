

void PileUpCorrection(TString infilename = "profiles_sys_min_pileup.root",TString histfile ="", TString outfilename = "PileUpCorrection.root")
{

  //  gSystem->Load("momentCode.so");
  gSystem->Load("femtoLibrary/momentCode.so");
  gSystem->Load("../lib/momentCode.so");

  PileUpCorrection * puCor = new PileUpCorrection();
  puCor->LoadCumulant(infilename,0);
  puCor->LoadPileUpHistograms(histfile);
  puCor->CorrectionForMultRange(20,300);

}



