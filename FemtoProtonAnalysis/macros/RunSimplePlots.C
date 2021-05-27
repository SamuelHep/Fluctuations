

void RunSimplePlots(TString filelist,TString parameterFile, TString InputFile_forSeed)
//void RunSimplePlots(TString filelist="../input_parameters/femto_nCharge_test.list",TString parameterFile="../input_parameters/parameter_n0p5_0_norm.list", TString InputFile_forSeed="/gpfs/mnt/gpfs01/star/pwg/sheppel/femtoDsts_2_15/652B7E50FD98B25D966180D8B6A7AC7E_0.fDst.root")
{

  gSystem->Load("momentCode.so");
  // gSystem->Load("../lib/momentCode.so");

  //Get Number from filename
  TObjArray * tokens = InputFile_forSeed.Tokenize("_");
  TString numName = tokens->At(tokens->GetEntries()-1)->GetName();
  numName.ReplaceAll("",".fDst.root");

  stringstream strNum(numName.Data());
  int seed = 0;
  strNum >> seed;
  cout << "Seed number = " << seed << endl;

  TChain * tc = GetTChainFromList(filelist,"fDst");
  long int nentries = tc->GetEntries();
   
  //read in the parameterlist
  InputParameterList pl = ReadInputFile(parameterFile);
  pl.PrintParameters();

  SimplePlots(
	      tc, //Event Chain
	      nentries, //Number of Events 
	      pl, // Parameter list
	      TString::Format("simplePlots_%i.root",seed) // name of outfile
	      ); //Only Calculates the moment profiles

}


