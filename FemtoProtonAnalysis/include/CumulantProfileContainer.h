#ifndef CUMULANT_PROF_CONTAIN_H
#define CUMULANT_PROF_CONTAIN_H

#include <vector>
#include <map>

#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"

class CumulantProfileContainer : public TObject {

 public:
  CumulantProfileContainer();
  CumulantProfileContainer(int iName);
  CumulantProfileContainer(TFile *f,int iName);
  CumulantProfileContainer(int iName,TString label,TH1F * multHist,LongDouble_t Cumulants[300][5],int MaxMult);

  ~CumulantProfileContainer();

  void FillProfile(int cent,std::vector<std::vector<double>> m_r_s);
  void FactorialCumulantsToNormalCumulants();
  void ReBinAllGraphs(std::vector<int> binEdges, std::vector<double> binLabels);

  void SetGraph(int n,TGraphErrors * gr){ _gr_vec[n] = (TGraphErrors*) gr->Clone(); }
  
  TGraphErrors * GetWeightGraph(){ return GetNGraph(0); }

  TGraphErrors * GetQ1Graph(){ return GetNGraph(1); }
  TGraphErrors * GetQ2Graph(){ return GetNGraph(2); }
  TGraphErrors * GetQ3Graph(){ return GetNGraph(3); }
  TGraphErrors * GetQ4Graph(){ return GetNGraph(4); }

  TGraphErrors * GetReBinnedQ1Graph(){ return GetNGraph(5); }
  TGraphErrors * GetReBinnedQ2Graph(){ return GetNGraph(6); }
  TGraphErrors * GetReBinnedQ3Graph(){ return GetNGraph(7); }
  TGraphErrors * GetReBinnedQ4Graph(){ return GetNGraph(8); }

  TGraphErrors * GetQ2_Q1_Graph(){ return GetNGraph(9); }
  TGraphErrors * GetQ3_Q2_Graph(){ return GetNGraph(10); }
  TGraphErrors * GetQ4_Q2_Graph(){ return GetNGraph(11); }

  TGraphErrors * GetReBinnedK1Graph(){ return GetNGraph(12); }
  TGraphErrors * GetReBinnedK2Graph(){ return GetNGraph(13); }
  TGraphErrors * GetReBinnedK3Graph(){ return GetNGraph(14); }
  TGraphErrors * GetReBinnedK4Graph(){ return GetNGraph(15); }

  TGraphErrors * GetK2_K1_Graph(){ return GetNGraph(16); }
  TGraphErrors * GetK3_K1_Graph(){ return GetNGraph(17); }
  TGraphErrors * GetK4_K1_Graph(){ return GetNGraph(18); }

  int NGraphs(){  if ( !_qGraphsFilled ) FactorialCumulantsToNormalCumulants(); return _gr_vec.size(); }
  TGraphErrors * GetNGraph(int n);
  
  std::vector<TProfile*> GetProfileVec(){ return _prof_vec; }
  std::vector<TGraphErrors*> GetGraphVec(){ return _gr_vec; }

  std::map<TString,double> InitMap(int cent);
  
 private:
  double Avg_M(int cent,int index);
  TString IndexToName(int i);
  void FillFactorialCumulantQuants(std::vector<std::vector<double>> & m_r_s);
  TGraphErrors * ReBin(std::vector<int> &binEdges,std::vector<double> &binLabels,TGraphErrors * Ngr, TGraphErrors * Qgr);
  TGraphErrors * GetQnGraph(int n);
  TGraphErrors * RatioGraphs(TGraphErrors * num_gr, TGraphErrors * denom_gr);
  
  TGraphErrors * CalcK(int i);
  TGraphErrors * CalcK1();
  TGraphErrors * CalcK2();
  TGraphErrors * CalcK3();
  TGraphErrors * CalcK4();

  //member vars
  TGraphErrors * _gr_N;
  TGraphErrors * _gr_Q1;
  TGraphErrors * _gr_Q2;
  TGraphErrors * _gr_Q3;
  TGraphErrors * _gr_Q4;
  TGraphErrors * _gr_binned_Q1;
  TGraphErrors * _gr_binned_Q2;
  TGraphErrors * _gr_binned_Q3;
  TGraphErrors * _gr_binned_Q4;
  TGraphErrors * _gr_Q2_Q1;
  TGraphErrors *   _gr_Q3_Q2;
  TGraphErrors *   _gr_Q4_Q2;

  TGraphErrors * _gr_binned_K1;
  TGraphErrors * _gr_binned_K2;
  TGraphErrors * _gr_binned_K3;
  TGraphErrors * _gr_binned_K4;
  TGraphErrors * _gr_K2_K1;
  TGraphErrors *   _gr_K3_K1;
  TGraphErrors *   _gr_K4_K1;


  std::vector<TGraphErrors*> _gr_vec;
  std::vector<TProfile*> _prof_vec;

  long double _m[23];
  bool _qGraphsFilled;
    
  ClassDef(CumulantProfileContainer,1)

  
};

#endif
