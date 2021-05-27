

TString name(int i, int j)
{

  TString outDir = "/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/out_sys_5_18/cumulants/";

  TString prefix = "out";
  TString l1[2] = {"n0p1_0","0_0p1"};
  //  TString l1[14] = {"n0p1_0p1","n0p3_0","n0p4_0","n0p5_0","n0p5_0_pt1",
  //		   "n0p5_0_pt2","n0p5_0_pt3","n0p7_0","n0p8_0","n0p9_0",
  //		   "n0p2_0","n0p3_n0p1","n0p4_n0p2","n0p5_n0p3"};
  TString l2[11] = {"norm","SYS1","SYS2","SYS3","SYS4","SYS5","SYS6","SYS7","SYS8","SYS9","SYS10"};

  if ( j==-1 ) return TString::Format("%scumulants_%s.root",outDir.Data(),l1[i].Data());
  else return TString::Format("%s%s_%s_%s.root",outDir.Data(),prefix.Data(),l1[i].Data(),l2[j].Data());

}


void RunSystematicCorrection()
{
  
  gSystem->Load("../lib/momentCode.so");

  cout << name(0,0) << endl;
  
  for ( int i=0; i<2;i++ )
    {
      SystematicCalculation * sys_calc = new SystematicCalculation();
      
      sys_calc->AddNominal( name(i,0) );
      sys_calc->AddSysPair( name(i,1), name(i,2) );
      sys_calc->AddSysPair( name(i,3), name(i,4) );
      sys_calc->AddSysPair( name(i,5), name(i,6) );
      sys_calc->AddSysPair( name(i,7), name(i,8) );
      //      sys_calc->AddSysPair( name(i,9), name(i,10) );
      
      sys_calc->Calculate();
      
      sys_calc->WriteToOutFile( name(i,-1) );
      
      delete sys_calc;
    }
  
}


