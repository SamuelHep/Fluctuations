

void MakeMultHists(TString filelist,TString parameterFile,char *index="0")
{

  //  gSystem->Load("../lib/momentCode.so");
  gSystem->Load("momentCode.so");

  //Get Number from filename
  TChain * tc = GetTChainFromList(filelist,"fDst");
  long int nentries = tc->GetEntries();
   
  //read in the parameterlist
  InputParameterList pl = ReadInputFile(parameterFile);
  pl.PrintParameters();

  TString outfilename = TString::Format("FxtMultHistograms_%s.root",index);
    
  MakeMultHists( 
		tc, //Event Chain
		nentries, //Number of Events 
		pl, // Parameter list
		outfilename // outfile name 
		 ); //Only Calculates the moment profiles

  

}



