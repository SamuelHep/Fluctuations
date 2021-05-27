#ifndef PROC_EMBEDDING_H
#define PROC_EMBEDDING_H

#include <map>
#include <vector>
#include "TString.h"
#include "InputParameterList.h"

class TTree;
class TH2D;
class TF1;
class TCanvas;
class TGraphErrors;
class TGraphAsymmErrors;

class ProcessEmbedding : public TObject {

 public:
  
  ProcessEmbedding(TString embed_filename,TString outfilename, TString particleName, InputParameterList & pl, TString sysName );
  ~ProcessEmbedding();

 private:
  Int_t Run();
  Int_t InitHistograms();
  Int_t RatioHistograms();
  TGraphAsymmErrors * RatioSlice(TString name,int ibin);
  Double_t Rapidity(double eta, double pT, double mass);

  void FitAllCentSlices();
  void PlotSlices();

  InputParameterList  _pl;
  TTree * tree;
  
  Float_t FXTMult3;
  Float_t vtx, vty, vtz;
  Float_t rcPt, rcP, rcEta, rcPhi, rcCharge, rcNfit, rcNposs;
  Float_t rcDca, rcDcaXy, rcDcaZ;
  Float_t mcPt, mcP, mcEta;

  Float_t particle_mass;

  TH2D * h_reco;
  TH2D * h_emb;
  TH2D * h_ratio;

  TString _sysName;

  std::map<TString,TH2D*> h2;
  std::map<TString,TH2D*> h2_ratio;
  std::map<TString,TGraphAsymmErrors*> h1;
  std::map<TString,TF1*> f1;
  std::map<TString,std::vector<TF1*>> slices;
  std::map<TString,TString> fit_name;

  std::vector<int> _bin_labels;
  std::vector<int> _bin_colors;
  std::vector<int> _bin_edges;
  std::vector<TCanvas*> _canvases;
  std::vector<TCanvas*> _canvases_slices;
  int _NCENTBINS;

  std::vector<TGraphErrors*> _gr_vec;
  std::vector<TGraphErrors*> _gr_par;
  std::vector<TF1*> _f_vec;
  TCanvas * _centCan;

  std::vector< TF1* > f1_vec;
  std::vector< TString > name_vec;


};


#endif
