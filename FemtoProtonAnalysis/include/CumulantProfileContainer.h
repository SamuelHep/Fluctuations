#ifndef CUMULANT_PROF_CONTAIN_H
#define CUMULANT_PROF_CONTAIN_H

#include <vector>
#include <map>

#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"

const static int kPrimary = -999;

class CumulantProfileContainer : public TObject {

 public:

  CumulantProfileContainer(int iName=-999);
  CumulantProfileContainer(int iName,TString label,TH1F * multHist,LongDouble_t Cumulants[300][7],int MaxMult);
  CumulantProfileContainer(TFile *f,int iName=-999);

  ~CumulantProfileContainer();

  void InitGraphs(int iName); //Initialize the _gr map
  void InitProfiles(int iName,TFile * file=NULL); //Initialize the _profile map 

  // Called by FillProfile() and generates the values to be stored in profiles 
  void FillBiVariateMoments(std::vector<std::vector<long double>> & m_r_s); 

  // Called once per event. Fills the profiles 
  void FillProfile(int cent,std::vector<std::vector<long double>> m_r_s);

  // Converts the bivariate moments to bivariate cumulants for a given centrality bin
  void InitCumulantMap(int cent); 

  // Converts the moment profiles into cumulant graphs 
  void MomentsToCumulants();
  
  // Calls the ReBin() method for every Cumulant, Factorial Cumulant and Ratio
  void ReBinAllGraphs(std::vector<int> binEdges, std::vector<double> binLabels);

  // Rebins a graph by bin edges and binlabels
  void ReBin(std::vector<int> &binEdges,std::vector<double> &binLabels,TGraphErrors * Ngr, TGraphErrors * Qgr,TGraphErrors * Qgr_cbwc);
  
  // Takes the ratios of two graphs 
  void RatioGraphs(TGraphErrors * num_gr, TGraphErrors * denom_gr, TGraphErrors * ratioGr);
  
  // Set a graph 
  void SetGraph(TString name,TGraphErrors * gr){ _gr[ name ] = (TGraphErrors*) gr->Clone(); }

  // Get Graphs in container
  TGraphErrors * GetWeightGraph(){ return _gr["N"]; }
  TGraphErrors * GetGraph(TString name){ return _gr[ name ]; }

  // Amend Graph Suffix.. adds the suffix to every graph
  void AmendGraphSuffix( TString suffix );

  // Number of graphs in graph map
  int NGraphs(){  _gr.size(); }
  
  std::map<TString,TProfile*> GetProfileMap(){ return _profile ; }
  std::map<TString,TGraphErrors*> GetGraphMap(){ return _gr    ; }
  
 private:
  double Avg_M(int cent,int index);
  TString IndexToName(int i);

  TGraphErrors * ReBin(std::vector<int> &binEdges,std::vector<double> &binLabels,TGraphErrors * Ngr, TGraphErrors * Qgr);
  TGraphErrors * GetQnGraph(int n);
  TGraphErrors * RatioGraphs(TGraphErrors * num_gr, TGraphErrors * denom_gr);
  
  std::map<TString,TGraphErrors*> _gr;
  std::map<TString,TProfile*> _profile;
  std::map<TString,long double> m;
  std::map<TString,long double> q;

  
    
  ClassDef(CumulantProfileContainer,1)

  
};

#endif
