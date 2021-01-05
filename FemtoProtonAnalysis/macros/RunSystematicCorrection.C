

TString name(int i, int j)
{

  TString outDir = "/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/cumulants/";

  TString prefix = "out";
  TString l1[14] = {"n0p1_0p1","n0p2_0p1","n0p3_0p1","n0p4_0p1","n0p5_0p1",
		   "n1_n0p5","n0p5_0","n0p5_0_pt1","n0p5_0_pt2","n0p5_0_pt3",
		   "n0p2_0","n0p3_n0p1","n0p4_n0p2","n0p5_n0p3"};
  TString l2[8] = {"norm","SYS1","SYS2","SYS3","SYS4","SYS5","SYS6","SYS7"};

  if ( j==-1 ) return TString::Format("%scumulants_%s.root",outDir.Data(),l1[i].Data());
  else return TString::Format("%s%s_%s_%s.root",outDir.Data(),prefix.Data(),l1[i].Data(),l2[j].Data());

}


void RunSystematicCorrection()
{
  
  gSystem->Load("../lib/momentCode.so");


  for ( int i=0; i<14;i++ )
    {
      SystematicCalculation * sys_calc = new SystematicCalculation();
      
      sys_calc->AddNominal( name(i,0) );
      //      sys_calc->AddSysSingle( name(i,1) );
      sys_calc->AddSysPair( name(i,2), name(i,3) );
      sys_calc->AddSysPair( name(i,4), name(i,5) );
      sys_calc->AddSysPair( name(i,6), name(i,7) );
      
      sys_calc->Calculate();
      
      sys_calc->WriteToOutFile( name(i,-1) );
      
      delete sys_calc;
    }

}


