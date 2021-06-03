


void RunProcessEmbedding(
			 TString input="/star/data03/pwg/bkimel/muEmbed/3GeV_Feb4/pro/muEmbed_Feb4_pro.root", 
			 TString analysisDir="",
			 TString particle="p"
			)
{
  
  gSystem->Load("../lib/momentCode.so");

  TString output[3];
  output[0] ="/rootfiles/eff/protonEmbedding_nhf10.root";
  output[1] ="/rootfiles/eff/protonEmbedding_nhf12.root";
  output[2] ="/rootfiles/eff/protonEmbedding_nhf15.root";


  TString params[3] = 
    {
      "/input_parameters/parameter_n0p5_0_norm.list",
      "/input_parameters/parameter_n0p5_0_SYS1.list",
      "/input_parameters/parameter_n0p5_0_SYS2.list"
    }

 for (int i=0;i<3;i++)
    {
      output[i]= TString::Format("%s%s",analysisDir.Data() , output[i].Data());
      params[i]= TString::Format("%s%s",analysisDir.Data() , params[i].Data());

      cout << output[i] << endl;
      cout << params[i] << endl;
    }

  InputParameterList pl[3];
 for (int i=0;i<3;i++)
    {
      pl[i] = ReadInputFile(params[i]);
    }

  ProcessEmbedding( input, output[0] , particle, pl[0],"nhf_10");
  ProcessEmbedding( input, output[1] , particle, pl[1],"nhf_12");
  ProcessEmbedding( input, output[2] , particle, pl[2],"nhf_15");

}
