#ifndef PROTON_EFFICIENCY_H
#define PROTON_EFFICIENCY_H

#include <vector>
#include "TObject.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TF1.h"
#include "InputParameterList.h"

#include <map>


class TGraph;
class TGraphAsymmErrors;

class ProtonEfficiency : public TObject
{

 public:
  ProtonEfficiency();
  ProtonEfficiency(TString filename1,TString filename2, InputParameterList &pl);
  ~ProtonEfficiency();

  double GetConstantEfficiency();
  void   SetConstantEfficiency(double val);
  double GetEff(int charge, double pt, double pz,bool tofmatch);

  std::map<TString,std::vector< TGraphAsymmErrors* >> _tg_tpc_eff; 
  std::map<TString,std::vector< TF1* >>  _tf1_tpc_eff; 
  std::map<TString,std::vector< double >>  _xMin; 

  std::map<TString,std::vector< TGraph* >> _tg_tof_eff; 

  std::vector<double> _tpc_bin_edges;
  std::vector<double> _tof_bin_edges;

  double Rapidity(double pt,double pz);
  double Pseudorapidity(double pt, double pz);

  double GetTPCEfficiency(double y, double pt);
  double GetTOFEfficiency(double y, double pt);

  void UnitTest(TH2D * tpcHist, TH2D * tofHist);
  
 private:

  int SetConstants();
  int LoadTPCEff(TString filename);
  int LoadTOFEff(TString filename);

  void SetLabels(InputParameterList &pl);

  bool constantEfficiencySet;
  double constantEfficiency;
  double _mass;
  
  std::vector<TString> tpcSysLabels;
  std::vector<TString> tofSysLabels;

  TString tofLabel, tpcLabel;

  TH2F * h_eff_tof;

  //  std::vector< TH1D* > _h1d_eff; 
  //  std::vector< TF1* >  _tf1_eff; 
  //  std::vector<double> _bin_edges;

  ClassDef(ProtonEfficiency,1)

};

#endif
