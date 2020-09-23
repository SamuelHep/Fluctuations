#ifndef PROTON_EFFICIENCY_H
#define PROTON_EFFICIENCY_H

#include <vector>
#include "TObject.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TF1.h"

class ProtonEfficiency : public TObject
{

 public:
  ProtonEfficiency();
  ProtonEfficiency(TString filename);
  ProtonEfficiency(TString filename1,TString filename2);
  ~ProtonEfficiency();

  double GetConstantEfficiency();
  void   SetConstantEfficiency(double val);
  double GetEff(int charge, double pt, double pz);
  double GetEff(int charge, double pt, double pz,bool tofmatch);

  std::vector< TH1D* > _h1d_eff; 
  std::vector< TF1* >  _tf1_eff; 
  std::vector<double> _bin_edges;

  double Rapidity(double pt,double pz);
  double GetTPCEfficiency(double y, double pt);
  double GetTOFEfficiency(double y, double pt);
  
 private:

  int SetConstants();
  int LoadTPCEff(TString filename);
  int LoadTOFEff(TString filename);

  bool constantEfficiencySet;
  double constantEfficiency;
  double _mass;
  
  TH2F * h_eff_tof;

  //  std::vector< TH1D* > _h1d_eff; 
  //  std::vector< TF1* >  _tf1_eff; 
  //  std::vector<double> _bin_edges;

  ClassDef(ProtonEfficiency,1)

};

#endif
