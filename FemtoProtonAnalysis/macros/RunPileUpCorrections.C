


void RunPileUpCorrections(TString profileDir, TString outDir, TString pileup_best, TString pileup_minus, TString pileup_plus)
{

  gSystem->Load("../lib/momentCode.so");
  TString histfilename[3] = { pileup_best, pileup_minus, pileup_plus };

  TString prefix="profiles_";
  TString out_prefix="out";
  TString label1[8] = {"n0p2_0","n0p3_0","n0p4_0","n0p5_0","n0p5_0_pt1","n0p5_0_pt2","n0p5_0_pt3","n0p5_0_pt4"};
  TString label2[9] = {"norm","SYS1","SYS2","SYS3","SYS4","SYS5","SYS6","SYS7","SYS8"};
  TString pileup_labels[2] = {"SYS9","SYS10"};
  TString cbwc_labels[2] = {"SYS11","SYS12"};

  for (int i=0;i<8;i++)
    {
      for (int j=0;j<9;j++)
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
