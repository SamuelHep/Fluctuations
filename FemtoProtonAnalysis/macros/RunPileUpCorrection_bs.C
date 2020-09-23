


//TString histfilename = "/star/u/yuzhang1/data01/3GeV_UrQMD/cas/Y-0.5_0_Y0.1/mix/Data/histos/histo_113.root",

void RunPileUpCorrection_bs(
			    TString infilename = "profiles_0w5_Effm5.root",
			    TString outfilename = "Corrected_0w5_Effm5_toyModel.root",
			    TString histfilename = "rootfiles/CorPlots_nocut.root",
			    Bool_t urqmdHists = 0
			    )
{

  gSystem->Load("../lib/momentCode.so");
  RunPileUpCorr(infilename,outfilename,histfilename,urqmdHists);

}
