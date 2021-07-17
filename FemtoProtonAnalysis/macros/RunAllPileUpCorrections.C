


//TString histfilename = "/star/u/yuzhang1/data01/3GeV_UrQMD/cas/Y-0.5_0_Y0.1/mix/Data/histos/histo_113.root",

void RunAllPileUpCorrections()
{

  gSystem->Load("../lib/momentCode.so");

  TString profileDir = "/star/u/sheppel/test_code/Fluctuations/FemtoProtonAnalysis/rootfiles/cumulants/";
  TString outDir = "/star/u/sheppel/test_code/Fluctuations/FemtoProtonAnalysis/rootfiles/cumulants/";
  TString pileupDir = "/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/pileupCor/";

  TString histfilename[3] = {"test_data_best.root","test_data_best_minus.root","test_data_best_plus.root"};

  for ( int i=0;i<3;i++)
    {
      histfilename[i] = pileupDir + histfilename[i];
    }

  TString prefix="profiles_";
  TString out_prefix="out";
  //  TString label1[4] = {"n0p2_0","n0p3_0","n0p4_0","n0p5_0"};
  //  TString label1[4] = {"n0p5_0_pt1","n0p5_0_pt2","n0p5_0_pt3","n0p5_0"};
  //  TString label1[4] = {"n0p1_0p1","n0p5_0_pt2","n0p5_0_pt3","n0p5_0"};
  //  TString label1[2] = {"n0p1_0","0_0p1"};
  TString label1[2] = {"","n0p5_0"};
  //  TString label2[10] = {"norm","SYS1","SYS2","SYS3","SYS4","SYS5","SYS6","SYS7","SYS8","SYS9"};
  TString label2[4] = {"","SYS_dca2p5","SYS_dca2","SYS_dca1"};
  //  TString pileup_labels[2] = {"SYS7","SYS8"};
  //  TString cbwc_labels[2] = {"SYS9","SYS10"};

  for (int i=1;i<2;i++)
    {
      for (int j=1;j<4;j++)
	{
	  TString infilename  = TString::Format("%s%s%s_%s.root",profileDir.Data(),prefix.Data(),label1[i].Data(),label2[j].Data());
	  TString outfilename = TString::Format("%s%s_%s_%s.root",outDir.Data(),out_prefix.Data(),label1[i].Data(),label2[j].Data());

	  cout << "INFILE =" << infilename << endl;
	  cout << "OUTFILE=" << outfilename << endl;
	  
	  RunPileUpCorr( infilename, outfilename, histfilename[0], 0);

	  if ( j==0 )
	    {
	      outfilename = TString::Format("%s%s_%s_%s.root",outDir.Data(),out_prefix.Data(),label1[i].Data(),pileup_labels[0].Data());
	      RunPileUpCorr( infilename, outfilename, histfilename[1], 0);
	      outfilename = TString::Format("%s%s_%s_%s.root",outDir.Data(),out_prefix.Data(),label1[i].Data(),pileup_labels[1].Data());
	      RunPileUpCorr( infilename, outfilename, histfilename[2], 0);

	      outfilename = TString::Format("%s%s_%s_%s.root",outDir.Data(),out_prefix.Data(),label1[i].Data(),cbwc_labels[0].Data());
	      RunPileUpCorr( infilename, outfilename, histfilename[0],0, -1);
	      outfilename = TString::Format("%s%s_%s_%s.root",outDir.Data(),out_prefix.Data(),label1[i].Data(),cbwc_labels[1].Data());
	      RunPileUpCorr( infilename, outfilename, histfilename[0],0, 1);
	    }

	}
    }

}
