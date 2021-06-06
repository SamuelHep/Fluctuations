
/*
This a crude bit of code to generate all of the parameter.list files for every pt, rapidity and sys configuration
Run root -q ParseFile.C
input_parameters folder will be populated with the correct input parameters
*/

void ParseFile()
{

  gSystem->Load("../lib/momentCode.so");

  ofstream base_input_file("../input_parameters/parameter_n0p2_0_norm.list");

  base_input_file << "Dca "              << 3 << endl;
  base_input_file << "EffMultiplier "    << 1 << endl;
  base_input_file << "Mass2High "        << 1.2 << endl;
  base_input_file << "Mass2Low "         << 0.6 << endl;
  base_input_file << "Mom "              << 2 << endl;
  base_input_file << "NHitsFitMin "      << 10 << endl;
  base_input_file << "NSigmaPionCut "    << 0 << endl;
  base_input_file << "NSigmaProtonCut " << 3 << endl;
  base_input_file << "PtHigh "           << 2 << endl;
  base_input_file << "PtLow "            << 0.4 << endl;
  base_input_file << "RapHigh "          << 1.049 << endl;
  base_input_file << "RapLow "           << 0.849 << endl;
  base_input_file << "VrMax "            << 1.5 << endl;
  base_input_file << "VzMax "            << 202 << endl;
  base_input_file << "VzMin "            << 199.5 << endl;


  ParseSysFile( "../input_parameters/parameter_n0p2_0_norm.list","../input_parameters/rap.list","n0p2_0_norm");
  ParseSysFile( "../input_parameters/parameter_n0p5_0_norm.list","../input_parameters/pt.list","norm");
    
  ParseSysFile( "../input_parameters/parameter_n0p2_0_norm.list","../input_parameters/sys.list","norm");
  ParseSysFile( "../input_parameters/parameter_n0p3_0_norm.list","../input_parameters/sys.list","norm");
  ParseSysFile( "../input_parameters/parameter_n0p4_0_norm.list","../input_parameters/sys.list","norm");
  ParseSysFile( "../input_parameters/parameter_n0p5_0_norm.list","../input_parameters/sys.list","norm");

  ParseSysFile( "../input_parameters/parameter_n0p5_0_pt1_norm.list","../input_parameters/sys.list","norm");
  ParseSysFile( "../input_parameters/parameter_n0p5_0_pt2_norm.list","../input_parameters/sys.list","norm");
  ParseSysFile( "../input_parameters/parameter_n0p5_0_pt3_norm.list","../input_parameters/sys.list","norm");
  ParseSysFile( "../input_parameters/parameter_n0p5_0_pt4_norm.list","../input_parameters/sys.list","norm");

}
