

void StudyTOF(TString infilename="/star/data01/pwg/sheppel/femtoAnalysis/tof_eff/out_eta/tofMatch_fineBinning.root",TString outfilename="TofSlices_Canvases.root")
{

    gSystem->Load("../lib/momentCode.so");
    TOFEfficiencyMaker * tof = new TOFEfficiencyMaker();
    tof->StudyTOFEfficiency(infilename,outfilename);


}
