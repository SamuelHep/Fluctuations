#ifndef _GLAUBER_MINIMIZER_
#define _GLAUBER_MINIMIZER_

#include <string>
#include <vector>
#include <iostream>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
//#include <TH1D.h>
//#include <TFile.h>
class TH1D;
class TFile;
class TGraph2D;
class GlauberClass;
class TString;
//class TString;

class GlauberMinimizer{
public:
  GlauberMinimizer(TH1D* a_measuredMultiplicity, TFile* a_glauberFile, double a_lowMultCut, double a_highMultCut);
  //~GlauberMinimizer();

  void FitGlauberModel(double a_initialMu, double a_initialK, double a_initialX = -1, int a_numEventOverride = -1, TString outdir="specifyanOutDir!",TString ofile="SPECIFYAFILENAME.root");

  //For Minimizer
  double chiSquaredFunct(const double* a_param);

  // normal / geometric / mixed
  void setChiSqrdMethod(TString a_chiSqrdMethod){m_chiSqrdMethod = a_chiSqrdMethod;};



private:


  TH1D* m_multiplicity;
  TH1D* m_simMult;
  GlauberClass* m_glauberEvent;
  int m_numGlauberEntries;
  vector<double> m_nPart;
  vector<double> m_nColl;
  double m_lowMultCut;
  double m_highMultCut;
  int m_lowMultCutBin;
  int m_highMultCutBin;
  TGraph2D* m_chiSqrdGraph;
  vector< TH1D* > m_simHistory;
  vector< vector< double > > m_simParamHistory; // Npp K X Chi^2

  TString m_chiSqrdMethod;

  ROOT::Math::Minimizer* m_minimizer;
  //ROOT::Math::Minuit2Minimizer* m_minimizer;
  

};





#endif
