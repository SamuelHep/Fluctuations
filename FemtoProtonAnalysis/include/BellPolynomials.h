#ifndef BELL_POLY_H
#define BELL_POLY_H

#include "TObject.h"
#include "TString.h"

class BellPolynomials : public TObject {

 public:

  BellPolynomials();
  ~BellPolynomials(){};
  
  B(int n, int x, int index);

 private:

  int bellVal1;
  int bellVal2;
  int bellVal3[3];
  int bellVal4[4][2];
  int bellVal5[5][2];
  int bellVal6[6][3];

  TString errorMessage;

  ClassDef( BellPolynomials, 1);

};
