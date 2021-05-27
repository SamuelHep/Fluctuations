#include <cstring>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <numeric>
#include <algorithm>
#include <sstream>
#include "TChain.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "InputParameterList.h"
#include "CumulantProfileContainer.h"
#include "analysisUtil.h"
#include "PileUpCorrection.h"

using namespace std;


float getProtonDedxMeanShift(float p)
{
  if ( p<0.2 || p>=2.6 ) return 0;

  float Mean[24] = {-0.56,-0.23,0.02,0.3,0.3,
		    0.33,0.53,0.48,0.52,0.55,
		    0.56,0.59,0.62,0.64,0.64,
		    0.64,0.64,0.64,0.65,0.68,
		    0.75,0.71,0.58,0.64};

  return Mean[int((p-0.2)*10.)];

}
  

int ComputeBootstrapErrors(CumulantProfileContainer * primaryCpc, vector<CumulantProfileContainer*> cpc_vec)
{

  for ( auto & key_graph_primary : primaryCpc->GetGraphMap() )
    {
      TString key = key_graph_primary.first;
      TGraphErrors * primary_graph = key_graph_primary.second;
      std::vector<TGraphErrors*> bootstrap_graphs;

      for ( auto & bootstrap_container : cpc_vec )
	{
	  bootstrap_graphs.push_back( (TGraphErrors*) bootstrap_container->GetGraph( key ) );
	}

      BootstrapGraph( primary_graph , bootstrap_graphs );

    }
  
  return 0;

}


int BootstrapGraph(TGraphErrors * gr, std::vector<TGraphErrors*> gr_vec)
{

  for (int i=0;i<gr->GetN();i++)
    {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];

      //      cout << "original (x,y) = " << x << ", " << y << endl;

      std::vector<double> msArray;
      
      for (auto & gr_bs : gr_vec)
	{

	  if (i >= gr->GetN()) continue;	 	  

	  double x_bs = gr_bs->GetX()[i];
	  double y_bs = gr_bs->GetY()[i];
	  //	  cout << "    bs (x , y) = " << x_bs << ", " << y_bs << endl;      
	  
	  if ( fabs(x - x_bs) > 0.001 ) continue;
	  msArray.push_back( pow((y - y_bs),2) );
	}

      double sum = accumulate(msArray.begin(),msArray.end(),0.0);
      int n = gr_vec.size();
      if ( n < 2 ) n = 2;
      double rms = sqrt( sum/( n-1 ));

      gr->SetPointError(i,0,rms);
      
    }

  return 0;
    
}


TChain * GetTChainFromList(TString listname, TString treename)
{

  TChain * chain = new TChain( treename.Data() );

  if (listname.Contains(".list") ||
      listname.Contains(".txt") ||
      listname.Contains(".input") )
    {
      ifstream fileList(listname);
      if (!fileList.is_open())
	{
	  fputs("ERROR OPENING FILELIST",stderr);
	  exit (EXIT_FAILURE);
	}

      string filename;
      while(fileList.good())
	{
	  getline(fileList,filename);
	  if ( filename.length() != 0 && filename.find(".root") != std::string::npos)
	    {
	      std::cout << "Adding " << filename << std::endl;
	      chain->Add( filename.c_str() );
	    }
	}
      fileList.close();
    }
  else
    {
      std::cout << "Adding " << listname.Data() << std::endl;
      chain->Add(listname.Data());
    }

  if (!chain)
    {
	  fputs("ERROR NO CHAIN",stderr);
	  exit (EXIT_FAILURE);
    }

  return chain;  
}


void ParseSysFile(TString nominal_file_name,TString list_of_sys_file_name, TString nominal_qualifier)
{

  ifstream list_of_sys( list_of_sys_file_name );
  if ( !list_of_sys )
    {
      cerr << "Unable to open sys file";
      exit(1);
    } 
  
  string line;
  while ( getline( list_of_sys, line ) )
    {
      istringstream iss(line);
      string syslabel, var;
      double value;
      if ( !(iss >> syslabel >> var >> value) ) {break;}
      GenerateSysFile( nominal_file_name , nominal_qualifier, syslabel, var, value );
    }

}

void GenerateSysFile(TString nominal_file_name, TString nominal_qualifier, TString sys_label, TString var, double value)
{

  //First get the nominal input parameter list
  InputParameterList input_par = ReadInputFile( nominal_file_name );
  
  //Apply Systematic cut
  if ( !input_par.ApplySystematic( var , value ) ) 
    {    
      cerr << "Unable to apply systematic to given variable";
      exit(1);
    }
      
  TString sys_name = nominal_file_name.ReplaceAll( nominal_qualifier , sys_label );

  WriteInputFile( input_par , sys_name );
 
}

void WriteInputFile(InputParameterList parList,TString outfilename)
{

  ofstream outFile(outfilename);
  for (auto &param : parList.GetParameterMap())
    {
      outFile << param.first << " " << param.second << endl;
    }
  outFile.close();
}

InputParameterList ReadInputFile(TString inputfile)
{

  ifstream inFile(inputfile);
  if (!inFile)
    {
      cerr << "Unable to open input parameter file";
      exit(1);
    }

  InputParameterList parList;
  
  string line;
  while (getline(inFile,line))
    {
      istringstream iss(line);
      string label;
      double val;
      if(!(iss >> label >> val)){ break; } //error
      parList.Read(label,val);
    }

  return parList;
}


void RunPileUpCorr(TString infilename, TString outfilename, TString histfilename, Bool_t urqmdHists,int shift_cut)
{
  
  cout << "INFILE  =" << infilename << endl;
  cout << "OUTFILE =" << outfilename << endl;
  cout << "HISTFILE=" << histfilename << endl;

  /*  std::vector<double> binLabels = { 10, 21, 32, 50, 74, 105, 149, 251};
  std::vector<int> binEdges     = {2, 5, 8, 13, 24, 29, 41, 50, 80};
  
  std::vector<double> binLabels = {39, 51, 64, 81, 100, 121, 146,  176, 210  , 248, 293, 345};
  std::vector<int> binEdges     = {2, 5, 8, 11,  13, 16, 20, 24, 29, 34, 41, 50, 90};
  */
  std::vector<double> binLabels = { 47, 70, 107, 157, 219, 282,  326};
  //CURRENT CENTRALITY DEFINITION
  std::vector<int> binEdges     = { 4 + shift_cut,    6 + shift_cut,    10 + shift_cut,    16 + shift_cut,   25 + shift_cut,  38 + shift_cut, 48 + shift_cut, 79 + shift_cut};

  //  std::vector<int> binEdges     = { 5,    7,    11,    18,   26,  40, 49,   90}; //Up one
  //  std::vector<int> binEdges     = { 3,    5,    9,    16,   24,  38,  47,   90}; //Down one

  std::vector<CumulantProfileContainer*> cpc_vec;
  std::vector<CumulantProfileContainer*> cpc_uncor_vec;

  cout << "Run over primary cumulant profiles" << endl;

  //Run the correction for primary
  PileUpCorrection * puCorr_prime = new PileUpCorrection();
  if (!urqmdHists) puCorr_prime->LoadPileUpHistograms( histfilename );
  else puCorr_prime->LoadURQMDHistograms( histfilename );
  CumulantProfileContainer * cpc_uncor_prime = puCorr_prime->LoadCumulant(infilename,kPrimary);
  CumulantProfileContainer * cpc_prime = puCorr_prime->CorrectionForMultRange(kPrimary,20,100);
  
  //  cpc_prime->SetGraph("N",cpc_uncor_prime->GetWeightGraph());
  cpc_prime->ReBinAllGraphs(binEdges,binLabels);  
  cpc_uncor_prime->ReBinAllGraphs(binEdges,binLabels);  

  cout << "Run over bootstrap cumulant profiles" << endl; 
  
  //Run the correction for the bootstraps
  for( int iBs=0;iBs<10;iBs++ )
    {
      PileUpCorrection * puCorr = new PileUpCorrection();
      if (!urqmdHists) puCorr->LoadPileUpHistograms( histfilename );
      else puCorr->LoadURQMDHistograms( histfilename );
      CumulantProfileContainer * cpc_uncor = puCorr->LoadCumulant(infilename,iBs);
      CumulantProfileContainer * cpc = puCorr->CorrectionForMultRange(iBs,20,100);
      
      cpc->SetGraph("N",cpc_uncor->GetWeightGraph());
      cpc->ReBinAllGraphs(binEdges,binLabels);  
      cpc_uncor->ReBinAllGraphs(binEdges,binLabels);  

      cpc_vec.push_back(cpc);
      cpc_uncor_vec.push_back(cpc_uncor);
    }

  ComputeBootstrapErrors(cpc_prime,cpc_vec);
  ComputeBootstrapErrors(cpc_uncor_prime,cpc_uncor_vec);

  cpc_uncor_prime->AmendGraphSuffix("_uncor");

  TFile * outfile = new TFile(outfilename,"recreate");

  for ( auto & key_gr : cpc_prime->GetGraphMap() )
    {
      key_gr.second->Write();
    }

  for ( auto & key_gr : cpc_uncor_prime->GetGraphMap() )
    {
      key_gr.second->Write();
    }

  delete cpc_prime;
  delete cpc_uncor_prime;

  outfile->Close();

}


void RunNoCorr(TString infilename, TString outfilename) 
{
  
  std::vector<double> binLabels = { 47, 70, 107, 157, 219, 282,  326};
  //CURRENT CENTRALITY DEFINITION
  std::vector<int> binEdges     = { 4 ,    6 ,    10 ,    16 ,   25 ,  38 , 48 , 79 };

  std::vector<CumulantProfileContainer*> cpc_uncor_vec;

  cout << "Run over primary cumulant profiles" << endl;

  TFile * f = new TFile(infilename,"read");

  //Run the correction for primary
  CumulantProfileContainer * cpc_uncor_prime = new CumulantProfileContainer(f,kPrimary);
  cpc_uncor_prime->MomentsToCumulants();

  cpc_uncor_prime->ReBinAllGraphs(binEdges,binLabels);  
  cout << "Run over bootstrap cumulant profiles" << endl;
  
  //Run the correction for the bootstraps
  for( int iBs=0;iBs<10;iBs++ )
    {
      CumulantProfileContainer * cpc_uncor = new CumulantProfileContainer(f,iBs);
      cpc_uncor->MomentsToCumulants();
      cpc_uncor->ReBinAllGraphs(binEdges,binLabels);  
      cpc_uncor_vec.push_back(cpc_uncor);

    }

  ComputeBootstrapErrors(cpc_uncor_prime,cpc_uncor_vec);
  cpc_uncor_prime->AmendGraphSuffix("_uncor");

  TFile * outfile = new TFile(outfilename,"recreate");

  for ( auto & key_gr : cpc_uncor_prime->GetGraphMap() )
    {
      key_gr.second->Write();
    }

  delete cpc_uncor_prime;

}

void RunNoCorrNoCBWC(TString infilename, TString outfilename) 
{
  
  //  std::vector<double> binLabels = { 47, 70, 107, 157, 219, 282,  326};
  std::vector<double> binLabels = { 326, 282, 219, 157, 107, 70,  47};
  //CURRENT CENTRALITY DEFINITION
  std::vector<int> binEdges     = { 0 ,    1 ,    2 ,   3 ,   4, 5 , 6 , 7 };

  std::vector<CumulantProfileContainer*> cpc_uncor_vec;

  cout << "Run over primary cumulant profiles" << endl;

  TFile * f = new TFile(infilename,"read");

  //Run the correction for primary
  CumulantProfileContainer * cpc_uncor_prime = new CumulantProfileContainer(f,kPrimary);
  cpc_uncor_prime->MomentsToCumulants();

  cpc_uncor_prime->ReBinAllGraphs(binEdges,binLabels);  
  cout << "Run over bootstrap cumulant profiles" << endl;
  
  //Run the correction for the bootstraps
  for( int iBs=0;iBs<10;iBs++ )
    {
      CumulantProfileContainer * cpc_uncor = new CumulantProfileContainer(f,iBs);
      cpc_uncor->MomentsToCumulants();
      cpc_uncor->ReBinAllGraphs(binEdges,binLabels);  
      cpc_uncor_vec.push_back(cpc_uncor);

    }

  ComputeBootstrapErrors(cpc_uncor_prime,cpc_uncor_vec);
  cpc_uncor_prime->AmendGraphSuffix("_uncor");

  TFile * outfile = new TFile(outfilename,"recreate");

  for ( auto & key_gr : cpc_uncor_prime->GetGraphMap() )
    {
      key_gr.second->Write();
    }

  delete cpc_uncor_prime;

}


pair<bool,bool> PileUpBadEvent(int fxt, double nmip,double pipdu)
{

  bool GoodEpd =true;
  bool GoodTof =true;

  int centBin = CentBin(fxt);
  if (centBin < 0 ) return make_pair(false,false);
  if (centBin > 9 ) centBin = 9;
  
  //  int epdCuts[11] = {73,62,59,55,51,46,41,35,30,25,16};
  //    int epdCuts[11] = {73,62,59,55,51,46,41,35,30,25,16};
  //  int epdCuts[8]    = {71,73,74,75,75,75,75,82};
  double epdCuts[10] = {97.75, 89.25, 87.75, 85.75, 83.25, 80.25, 77.25, 74.25, 71.75, 68.75};
  double tofCuts[10] = {29.4098, 24.3502, 18.8585, 14.0175, 11.1732, 8.18495, 6.82341, 5.88937, 4.19889,0};
  if (nmip > epdCuts[centBin]) GoodEpd = false;

  //  if (centBin > 8 ) centBin = 8;
  if (pipdu < tofCuts[centBin]) GoodTof = false;

  return make_pair(GoodEpd,GoodTof);  
}

pair<bool,bool> PileUpBadEvent(int gmult, int fxt,double nmip,double pipdu)
{

  bool GoodEpd =true;
  bool GoodTof =false;

  int centBin = CentBinTofMatch(gmult);
  if (centBin < 0 ) return make_pair(false,false);
  if (centBin > 9 ) centBin = 9;
  
  double epdCuts[10] = {97.75, 89.25, 87.75, 85.75, 83.25, 80.25, 77.25, 74.25, 71.75, 68.75};
  if (nmip > epdCuts[centBin]) GoodEpd = false;

  double p0 = -4.56;
  double p1 = 0.0504;
  double p2 =  0.000816;

  double cut = p0 + p1*(fxt) + p2*(fxt*fxt);

  if (pipdu > cut && pipdu > 0) GoodTof = true;

  return make_pair(GoodEpd,GoodTof);  
}

pair<bool,bool> PileUpBadEventVariable(int fxt,double nmip,double pipdu,int index)
{
  
  int pileup_cut = 200;

  if ( index == 0 ) pileup_cut = 175; 
  if ( index == 1 ) pileup_cut = 180; 
  if ( index == 2 ) pileup_cut = 185; 
  if ( index == 3 ) pileup_cut = 190; 
  if ( index == 4 ) pileup_cut = 195; 
  if ( index == 5 ) pileup_cut = 200; 
  if ( index == 6 ) pileup_cut = 205; 
  if ( index == 7 ) pileup_cut = 210; 
  if ( index == 8 ) pileup_cut = 215; 
  if ( index == 9 ) pileup_cut = 220; 

  if ( fxt > pileup_cut ) return make_pair(false,false);
  return make_pair(true,true);

  /*

  bool GoodEpd =true;
  bool GoodTof =false;

  int centBin = CentBin(fxt);
  if (centBin < 0 && index != 9) return make_pair(false,false);
  
  if ( index >= 8  ) return make_pair(true,true); //if index is 0 or greater than 8, just say good event

  double epdCuts1[16]       = {41, 42, 43, 45, 49, 52, 56, 61, 65, 68, 71, 74, 76, 77, 80, 80};
  double epdCuts2[16]       = {51, 51, 52, 54, 58, 62, 66, 69, 73, 76, 80, 83, 85, 86, 90, 90};
  double epdCuts3[16]       = {74, 75, 76, 77, 80, 83, 86, 89, 92, 95, 98, 101, 104, 106, 110, 110};

  if (index == 0 || index == 2 || index == 3 ) //For first 3, use 80 cut
    {
      if ( nmip > epdCuts2[14 - centBin] ) GoodEpd = false;
    }
  if ( index == 1 ) //For 4, use 90 cut
    {
      if ( nmip > epdCuts1[14 - centBin] ) GoodEpd = false;
    }
  if ( index == 4 ) //For 5, use 50 cut
    {
      if ( nmip > epdCuts3[14 - centBin] ) GoodEpd = false;
    }

  double p0[3] = {-2.247,-3.03665,-3.82582};
  double p1[3] = {0.08050,0.02068,-0.03919};
  double p2[3] = {0.001,0.00117347,0.001344};

  int tIndex = 0;

  if ( index == 1 || index == 2 || index == 4 || index == 6 ) tIndex = 1;
  if ( index == 0 || index == 5 ) tIndex = 0;
  if ( index == 3 || index == 7 ) tIndex = 2;

  double cut = p0[tIndex] + p1[tIndex]*(fxt) + p2[tIndex]*(fxt*fxt);

  if (pipdu > cut && pipdu > 0) GoodTof = true;

  return make_pair(GoodEpd,GoodTof);  

  */

}


int CentBin(int mult)
{

  if (mult > 200) return -1;
  if (mult > 136) return 0;
  if (mult > 114) return 1;
  if (mult > 94) return 2;
  if (mult > 78) return 3;
  if (mult > 64) return 4;
  if (mult > 52) return 5;
  if (mult > 44) return 6;
  if (mult > 34) return 7;
  if (mult > 28) return 8;
  if (mult > 22) return 9;
  if (mult > 16) return 10;
  if (mult > 12) return 11;
  if (mult > 10) return 12;
  if (mult > 6) return 13;
  else return 14;
}

int CentBinTofMatch(int mult)
{

  if (mult > 80) return -1;
  if (mult > 48) return 0;
  if (mult > 40) return 1;
  if (mult > 34) return 2;
  if (mult > 28) return 3;
  if (mult > 22) return 4;
  if (mult > 18) return 5;
  if (mult > 14) return 6;
  if (mult > 12) return 7;
  if (mult > 10) return 8;
  if (mult > 8) return 9;
  if (mult > 6) return 10;
  if (mult > 4) return 11;
  else return 12;
}


int CentBin3(int mult)
{
  if (mult > 79) return -1;
  if (mult > 48) return 0;
  if (mult > 38) return 1;
  if (mult > 25) return 2;
  if (mult > 16) return 3;
  if (mult > 10) return 4;
  if (mult > 6) return 5;
  if (mult > 4) return 6;
  if (mult >= 0) return 7;
  //  if (mult > 14) return 7;
  //  if (mult > 11) return 8;
  //  if (mult > 9) return 9;
  //  if (mult > 7) return 10;
  //  if (mult > 6) return 11;
  //  if (mult > 4) return 12;
  //  if (mult > 3) return 13;
  else return -1;
}


int GetRunIndex( std::vector<int> &runvec, int runNum)
{

  auto it = std::find(runvec.begin(),runvec.end(),runNum);
  if (it != runvec.end())
    {
      return std::distance(runvec.begin(),it);
    }
  else
    {
      runvec.push_back(runNum);
      return std::distance(runvec.begin(),runvec.end()) - 1;
    }
}
