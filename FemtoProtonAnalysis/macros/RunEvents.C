

void RunEvents(TString filelist,TString parameterFile, TString InputFile_forSeed, TString tpc_efficiency_file, TString tof_efficiency_file)
{

  gSystem->Load("momentCode.so");
  //gSystem->Load("../lib/momentCode.so");

  //Get Number from filename
  TObjArray * tokens = InputFile_forSeed.Tokenize("_");
  TString numName = tokens->At(tokens->GetEntries()-1)->GetName();
  numName.ReplaceAll(".fDst.root","");

  cout << numName << endl;

  stringstream strNum(numName.Data());
  int seed = 0;
  strNum >> seed;
  cout << "Seed number = " << seed << endl;

  TChain * tc = GetTChainFromList(filelist,"fDst");
  long int nentries = tc->GetEntries();
   
  //read in the parameterlist
  InputParameterList pl = ReadInputFile(parameterFile);
  pl.PrintParameters();

  TString filename = parameterFile;
  filename.ReplaceAll("parameter_","");
  filename.ReplaceAll(".list","");
  filename = TString::Format("%s_profiles_%i.root",filename.Data(),seed); // name of outfile
  cout << "Filename = " << filename << endl;

  //Make Proton Efficiency Object 
  ProtonEfficiency * eff =  new ProtonEfficiency(tpc_efficiency_file, tof_efficiency_file, pl);
  //  ProtonEfficiency * eff =  new ProtonEfficiency("tpc_efficiency.root","tof_efficiency_fineBinning.root",pl);
  //  ProtonEfficiency * eff =  new ProtonEfficiency();
  //  eff->SetConstantEfficiency(1.0);
  
  //Make a trandom3 to pass random numbers to eventloop
  TRandom3 * rand = new TRandom3(seed);
    
  MomentEventLoopPrintToFile( 
			     tc, //Event Chain
			     nentries, //Number of Events 
			     pl, // Parameter list
			     eff,  // efficiency 
			     rand, // trandom3
			     10, // number of bootstraps
			     filename // outfile name 
			      ); //Only Calculates the moment profiles

}



