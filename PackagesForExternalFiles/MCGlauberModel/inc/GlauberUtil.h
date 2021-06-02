#ifndef GLAUBER_UTIL_H
#define GLAUBER_UTIL_H

#include "GlauberClass.h"

int HISTCOUNTER=0;

TH1D * MakeNegativeBinomialHist(Double_t npp,Double_t k,Double_t hardness,Int_t index,GlauberClass& glauberEvent);

#endif
