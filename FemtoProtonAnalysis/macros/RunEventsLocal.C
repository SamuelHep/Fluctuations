

void RunEventsLocal(TString filelist,TString parameterFile, TString outfilename)
{

  //  gSystem->Load("momentCode.so");
  gSystem->Load("../lib/momentCode.so");

  TChain * tc = GetTChainFromList(filelist,"fDst");
  long int nentries = tc->GetEntries();
   
  //read in the parameterlist
  InputParameterList pl = ReadInputFile(parameterFile);
  pl.PrintParameters();

  //Make Proton Efficiency Object 
  //ProtonEfficiency * eff =  new ProtonEfficiency("fitEmbedding.root","tof_eff.root");
  ProtonEfficiency * eff =  new ProtonEfficiency();
  eff->SetConstantEfficiency(1.0);
  
  //Make a trandom3 to pass random numbers to eventloop
  TRandom3 * rand = new TRandom3(123);
    
  MomentEventLoopLocal( 
		       tc, //Event Chain
		       nentries, //Number of Events 
		       pl, // Parameter list
		       eff,  // efficiency 
		       rand, // trandom3
		       10, // number of bootstraps
		       outfilename // outfile name 
			); //Only Calculates the moment profiles

  

}



