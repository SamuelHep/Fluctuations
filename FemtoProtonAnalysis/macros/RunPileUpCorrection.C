


void RunPileUpCorrection(TString infilename = "profiles_n0p5_0_norm.root",TString outfilename = "weighted_correction_n0p5_0_norm.root")
{
  
  gSystem->Load("../lib/momentCode.so");

  TString profileDir = "/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/profiles/newCode/";
  TString outDir = "/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/cumulants/newCode/";
  TString pileupDir = "/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/pileupCor/";

  TString histfilename[3] = {"test_data_best.root","test_data_best_minus.root","test_data_best_plus.root"};

  TString infile_path  = TString::Format("%s%s",profileDir.Data(), infilename.Data() );
  TString outfile_path = TString::Format("%s%s",outDir.Data(), outfilename.Data() );
  
  TString histcorr_path = TString::Format("%s%s", pileupDir.Data(), histfilename[0].Data() );

  RunPileUpCorr(infile_path, outfile_path, histcorr_path, 0);


}
