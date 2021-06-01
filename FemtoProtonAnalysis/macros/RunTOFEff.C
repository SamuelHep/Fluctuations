


void RunTOFEff(TString filelist, TString parameterFile,char *index="0")
{

  gSystem->Load("momentCode.so");

  // Make TChain
  TChain * tc = GetTChainFromList(filelist,"fDst");
  long int nentries = tc->GetEntries();

  //read in the parameter list
  InputParameterList pl = ReadInputFile(parameterFile);
  pl.PrintParameters();
 
  TString outfilename = TString::Format("tofMatch_%s.root",index);

  TOFEfficiencyMaker * tofMaker = new TOFEfficiencyMaker();
  tofMaker->LoopEvents(
		       tc, //Event Chain
		       nentries, //Number of Events 
		       pl, // Parameter list
		       outfilename // outfile name 
		       ); //Only Calculates the moment profiles
  
}
