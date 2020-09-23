R__LOAD_LIBRARY(momentCodeLibrary/momentCode.so);

InputParameterList ReadInputFile(TString inputfile);

void RunCumulant(TString inputfile="parameter.list",TString filelist="run_fxt.list")
{

  clock_t t = clock();
  
  //init the event chain and set the event branch address
  TChain * tc = GetTChainFromList(filelist,"fDst");
  //  long int nentries = tc->GetEntries();
  long int nentries = 10000000;
  
  //read in the parameterlist
  InputParameterList pl = ReadInputFile(inputfile);

  pl.PrintParameters();
  //Make CumulantContainer Object to store result
  CumulantProfileContainer * cpc = new CumulantProfileContainer();

  //Make Proton Efficiency Object 
  ProtonEfficiency * eff =  new ProtonEfficiency();
  //ProtonEfficiency * eff =  new ProtonEfficiency("../Eff_Proton_embedding_2D.root","../QA/QA_outfile_eff.root");
  eff->SetConstantEfficiency(1.0);

  //  cout << "test " << eff->GetEff(1,0.5,2.0) << endl;
  
  //Make a trandom3 to pass random numbers to eventloop
  TRandom3 * rand = new TRandom3();
    
  //  EventLoop( tc , nentries , pl , cpc , eff );

  vector<TGraphErrors*> gr; // Graph for Q1 - Q4
  vector<TGraphErrors*> gr1; // Graph for Q1 - Q4
  vector<TGraphErrors*> gr2; // Graph for Q1 - Q4
  vector<TProfile*> prof_vec;

  //Define the binning
  std::vector<double> binLabels = {118,170,240,308 ,395};
  std::vector<int> binEdges =  {12,19,27,39,47 ,80};  
  std::vector< CumulantProfileContainer* > cpc_vec;

  for (int i=0;i<20;i++)
    {
      cpc_vec.push_back(new CumulantProfileContainer(i));
    }
  
  EventLoopBootstrap( tc , nentries , pl , cpc, cpc_vec , eff , rand);

  cpc->ReBinAllGraphs(binEdges,binLabels);
  
  for (auto &v : cpc_vec)
    {
      v->ReBinAllGraphs(binEdges,binLabels);
    }

  ComputeBootstrapErrors(cpc,cpc_vec);

  prof_vec = cpc->GetProfileVec();
  gr = cpc->GetGraphVec();
  gr1 = cpc_vec[0]->GetGraphVec();
  gr2 = cpc_vec[1]->GetGraphVec();
  
  
  //for ( int i=0; i< cpc->GetNGraphs();i++)
  //    {
  //      gr.push_back( cpc->GetNGraph(i) );
  //    }
  
  TFile * outfile = new TFile("cumulants_outfile_epdtofcuts_dcacut_mom2GeV_moreBootstraps.root","recreate");
  outfile->cd();

  for (auto &g : gr)
    {
      g->Write();
    }

  for (auto &g : gr1)
    {
      g->Write();
    }

    for (auto &g : gr2)
    {
      g->Write();
    }

  
  for (auto &p :prof_vec)
    {
      p->Write();
    }

  t = clock() - t;
  cout << "Run Time = " << ((float)t/CLOCKS_PER_SEC) << " secs" << endl;


}


InputParameterList ReadInputFile(TString inputfile)
{

  ifstream inFile(inputfile);
  if (!inFile)
    {
      cerr << "Unable to open input parameter file";
      exit(1);
    }

  InputParameterList parList;
  
  string line;
  while (getline(inFile,line))
    {
      istringstream iss(line);
      string label;
      double val;
      if(!(iss >> label >> val)){ break; } //error
      parList.Read(label,val);
    }

  return parList;
}
