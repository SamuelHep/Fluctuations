

TString name(TString outDir, int i, int j)
{
  TString prefix = "out";
  TString l1[8] = {"n0p2_0","n0p3_0","n0p4_0","n0p5_0",
		    "n0p5_0_pt1","n0p5_0_pt2","n0p5_0_pt3","n0p5_0_pt4"};

  TString l2[11] = {"norm","SYS1","SYS2","SYS3","SYS4","SYS5","SYS6","SYS7","SYS8","SYS9","SYS10"};

  if ( j==-1 ) return TString::Format("%scumulants_%s.root",outDir.Data(),l1[i].Data());
  else return TString::Format("%s%s_%s_%s.root",outDir.Data(),prefix.Data(),l1[i].Data(),l2[j].Data());
}


void CalcSystematic(TString cumu_dir)
{
  
  gSystem->Load("../lib/momentCode.so");

  cout << name(cumu_dir,0,0) << endl;
  
  for ( int i=0; i<8;i++ )
    {
      SystematicCalculation * sys_calc = new SystematicCalculation();
      
      sys_calc->AddNominal( name(cumu_dir,i,0) );
      sys_calc->AddSysPair( name(cumu_dir,i,1), name(cumu_dir,i,2) );
      sys_calc->AddSysPair( name(cumu_dir,i,3), name(cumu_dir,i,4) );
      sys_calc->AddSysPair( name(cumu_dir,i,5), name(cumu_dir,i,6) );
      sys_calc->AddSysPair( name(cumu_dir,i,7), name(cumu_dir,i,8) );
      
      sys_calc->Calculate();
      
      sys_calc->WriteToOutFile( name(cumu_dir,i,-1) );
      cout << "done!" << endl;
      cout << endl;
      
      delete sys_calc;
    }
  
}


