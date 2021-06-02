
#include "TH1D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TNtuple.h"
#include <vector>
#include <iostream>

#include "ReturnStructs.h"

using namespace std;

//___________________________________________
void Delete(NegBinomialSearchResults *result){

  //This Function is responsible for deleting
  //the object passed into it
  delete result->bestFitHisto;
  delete result->chi2Graph1D;
  delete result->chi2Graph2D;
  delete result->Prob1D;
  delete result->nTup3D;
}

//____________________________________________
ClassImp(RefMultCorrInfo);

RefMultCorrInfo::RefMultCorrInfo(){

}

RefMultCorrInfo::~RefMultCorrInfo(){

}
