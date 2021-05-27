#ifndef TOF_EFF_MAKER_H
#define TOF_EFF_MAKER_H

#include "TGraphAsymmErrors.h"
#include "TChain.h"
#include "InputParameterList.h"
#include "TString.h"
#include <map>
#include "TF1.h"
#include "TCanvas.h"

class StFemtoEvent;
class TH3D;
class TGraph;
class TGraphAssymErrors;


class TOFEfficiencyMaker : public TObject {

 public:
  
  TOFEfficiencyMaker();
  ~TOFEfficiencyMaker(){};
  
  //Use to generate Efficiency plots, looping over fDsts
  void LoopEvents(TChain * tc,
		  long int entries,
		  InputParameterList & pl,
		  TString outfilename
		  );

  void StudyTOFEfficiency(TString filename,TString outfilename);
  
 private:
  
  StFemtoEvent * event;
  std::map<TString,TH3D*> TOFHists;

  std::vector<TCanvas*> PlotSlices(std::vector<TGraph*> &f, std::vector<TGraphAsymmErrors*> &gr, TString canName);
  
  double const mass;

  int const MAXFXT3;
  int const PTBINS;
  int const ETABINS;

  std::vector<TString> hist_names;



};

#endif
