
void RunSystematicEvents(TString filelist,TString parameterFile, TString InputFile_forSeed)
{

  gSystem->Load("momentCode.so");
  //gSystem->Load("femtoLibrary/momentCode.so");

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

  //Make Proton Efficiency Object 
  ProtonEfficiency * eff =  new ProtonEfficiency("fitEmbedding.root","tof_eff.root");
  //  eff->SetConstantEfficiency(1.0);
  
  //Make a trandom3 to pass random numbers to eventloop
  TRandom3 * rand = new TRandom3(seed);
    
  SystematicEventLoopPrintToFile( 
				 tc, //Event Chain
				 nentries, //Number of Events 
				 pl, // Parameter list
				 eff,  // efficiency 
				 rand, // trandom3
				 20, // number of bootstraps
				 "rap",
				 TString::Format("profiles_sys_%i.root",seed) // name of outfile
				  ); //Only Calculates the moment profiles

}



