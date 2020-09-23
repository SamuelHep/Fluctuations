


void ComputeSysAndErrors(TString infilename = "rapScan_profiles.root", TString outfilename = "check_rapScan.root")
{

  //    gSystem->Load("momentCode.so");
  gSystem->Load("../lib/momentCode.so");

    GetProfilesAndCalculateSystematics(infilename,outfilename,20,"rap",1.0,0.5);
  //  GetProfilesAndCalculateSystematics(infilename,outfilename,20,"mom",0.5,0.2);

}



