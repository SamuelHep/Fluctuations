


void RunPileUpCorrection(TString infilename = "profiles_noPUCut.root",TString outfilename = "test_correction.root")
{
  
  gSystem->Load("femtoLibrary/momentCode.so");

  //Bin labels to apply
  std::vector<double> binLabels;
  double binLabelArray[5] = {118,170,240,308 ,395};
  std::vector<int> binEdges;
  int binEdgeArray[6] = {12,19,27,39,47 ,80};
  for (int i=0;i<5;i++) binLabels.push_back(binLabelArray[i]);
  for (int i=0;i<6;i++) binEdges.push_back(binEdgeArray[i]);

  PileUpCorrection * puCorr = new PileUpCorrection();
  CumulantProfileContainer * cpc_uncor = puCorr->LoadCumulant(infilename,1);
  //puCorr->LoadPileUpHistograms( "FxtMult3_hists.root" );
  puCorr->LoadURQMDHistograms( "yHisto0.root" );
  CumulantProfileContainer * cpc = puCorr->CorrectionForMultRange(1,0,90);

  cpc->SetGraph(0,cpc_uncor->GetWeightGraph());
  cpc->ReBinAllGraphs(binEdges,binLabels);  
  cpc_uncor->ReBinAllGraphs(binEdges,binLabels);  

  TFile * outfile = new TFile(outfilename,"recreate");

  for( int i=0;i<cpc->NGraphs();i++)
    {
      TGraphErrors * gr = cpc->GetNGraph(i);
      TGraphErrors * gr_uncor = cpc_uncor->GetNGraph(i);
      gr->Write();
      gr_uncor->Write();
      delete gr;
      delete gr_uncor;
    }

}
