#ifndef SYS_CALC_H
#define SYS_CALC_H

#include "TObject.h"
#include "TString.h"

#include <vector>
#include "TGraphErrors.h"
#include <utility>

class SystematicCalculation : public TObject  {

 private:

  std::vector<TGraphErrors*>  _nominal_gr_vec;
  std::vector<TGraphErrors*>  _sys_gr_vec;
  std::vector<std::vector<TGraphErrors*>>  _nsys_single_gr;
  std::vector<std::vector<TGraphErrors*>>  _nsys_double_gr;
  std::vector<TString> _simple_names_norm;
  std::vector<TString> _simple_names_sys;

  std::vector<TString> _grnames;

 public:
  
  SystematicCalculation();
  ~SystematicCalculation();

  bool AddNominal(TString filename);
  bool AddSysSingle(TString filename);
  bool AddSysPair(TString filename1,TString filename2);
  bool Calculate();
  bool CheckGraphVectors();

  bool WriteToOutFile(TString outfilename,bool rename=true);

  ClassDef( SystematicCalculation, 1 );

};



#endif
