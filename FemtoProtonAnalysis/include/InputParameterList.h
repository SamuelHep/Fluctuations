#ifndef IMPUT_PARAMETER_H
#define IMPUT_PARAMETER_H

#include "TObject.h"
#include "TString.h"
#include <map>

class InputParameterList : public TObject
{

 public:

  InputParameterList();
  ~InputParameterList();

  void Read(std::string label,double val);
  void PrintParameters();
  double PtLow();
  double PtHigh();
  double RapidityLow();
  double RapidityHigh();
  double NHitsFitMin();
  double VzMin();
  double VzMax();
  double VrMax();
  double NSigmaProtonCut();
  double NSigmaPionCut();
  double Mass2High();
  double Mass2Low();
  double Mom();
  double EffMultiplier();
  double Dca();
  
 private:

  std::map<TString,double> _parameterMap;
  
  ClassDef(InputParameterList, 1)

    };

#endif
