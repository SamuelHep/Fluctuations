#include <map>
#include <iostream>
#include "InputParameterList.h"

using namespace std;

ClassImp(InputParameterList)
  
//_______________________________________
InputParameterList::InputParameterList()
{

  _parameterMap["PtLow"] = 0.4;
  _parameterMap["PtHigh"] = 2.0;
  _parameterMap["RapLow"] = -0.5;
  _parameterMap["RapHigh"] = 0.5;
  _parameterMap["NHitsFitMin"] = 10;
  _parameterMap["VzMin"] = 198;
  _parameterMap["VzMax"] = 202;
  _parameterMap["VrMax"] = 1.5;
  _parameterMap["NSigmaProtonCut"] = 2;
  _parameterMap["NSigmaPionCut"] = 0;
  _parameterMap["Mass2High"] = -1;
  _parameterMap["Mass2Low"] = -1;
  _parameterMap["Mom"] = 5;
  _parameterMap["EffMultiplier"] = 1;
  _parameterMap["Dca"] = 3.0;
  
}

//_______________________________________
InputParameterList::~InputParameterList()
{
}

void InputParameterList::PrintParameters()
{

  for (auto &v : _parameterMap)
    {
      cout << " " << v.first << " = " << v.second << endl;
    }

}

//_______________________________________
void InputParameterList::Read(std::string label,double val)
{
  TString str = static_cast<TString>(label);
  _parameterMap[str] = val; 
}

//_______________________________________
bool InputParameterList::ApplySystematic(TString var, double value)
{

  if ( !var.CompareTo("NHitsFitMin") )
    {
      _parameterMap["NHitsFitMin"] = value;
      return true;
    }

  if ( !var.CompareTo("EffMultiplier") )
    {
      _parameterMap["EffMultiplier"] = value;
      return true;
    }

  if ( !var.CompareTo("Mass2Shift") )
    {
      _parameterMap["Mass2High"] = _parameterMap["Mass2High"] + value;
      _parameterMap["Mass2Low"]  = _parameterMap["Mass2Low"] + value;
      return true;
    }

  if ( !var.CompareTo("RapLowShift") )
    {
      _parameterMap["RapLow"] = _parameterMap["RapLow"] + value;
      return true;
    }

  if ( !var.CompareTo("RapShift") )
    {
      _parameterMap["RapLow"] = _parameterMap["RapLow"] + value;
      _parameterMap["RapHigh"] = _parameterMap["RapHigh"] + value;
      return true;
    }

  if ( !var.CompareTo("PtHighShift") )
    {
      _parameterMap["PtHigh"] = _parameterMap["PtHigh"] + value;
      return true;
    }

  return false;

}

//_______________________________________
double InputParameterList::PtLow()
{
  return _parameterMap["PtLow"];
}

//_______________________________________
double InputParameterList::PtHigh()
{
  return _parameterMap["PtHigh"];
}

//_______________________________________
double InputParameterList::RapidityLow()
{
  return _parameterMap["RapLow"];
}

//_______________________________________
double InputParameterList::RapidityHigh()
{
  return _parameterMap["RapHigh"];
}

//_______________________________________
double InputParameterList::NHitsFitMin()
{
  return _parameterMap["NHitsFitMin"];
}

//_______________________________________
double InputParameterList::VzMin()
{
  return _parameterMap["VzMin"];
}

//_______________________________________
double InputParameterList::VzMax()
{
  return _parameterMap["VzMax"];
}

//_______________________________________
double InputParameterList::VrMax()
{
  return _parameterMap["VrMax"];
}

//_______________________________________
double InputParameterList::NSigmaProtonCut()
{
  return _parameterMap["NSigmaProtonCut"];
}

//_______________________________________
double InputParameterList::NSigmaPionCut()
{
  return _parameterMap["NSigmaPionCut"];
}

//_______________________________________
double InputParameterList::Mass2High()
{
  return _parameterMap["Mass2High"];
}

//_______________________________________
double InputParameterList::Mass2Low()
{
  return _parameterMap["Mass2Low"];
}

//_______________________________________
double InputParameterList::Mom()
{
  return _parameterMap["Mom"];
}

//_______________________________________
double InputParameterList::EffMultiplier()
{
  return _parameterMap["EffMultiplier"];
}

//_______________________________________
double InputParameterList::Dca()
{
  return _parameterMap["Dca"];
}





