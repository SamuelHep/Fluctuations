#include "TH1D.h"
#include "TString.h"
#include "TF1.h"
#include "TRandom3.h"
#include "GlauberClass.h"
#include "GlauberUtil.h"


TH1D * MakeNegativeBinomialHist(Double_t npp,Double_t k,Double_t hardness,Int_t index,GlauberClass& glauberEvent)
{

  TString name = TString::Format("negativeBinomialHist_%i",++HISTCOUNTER);
  TH1D * nbdHist = new TH1D(name,"",100,0,100); 
  TF1 * nbdDist = glauberEvent.GetNegativeBinomialDistribution();

  for ( int i=0;i<100;i++)
    {
      nbdHist->SetBinContent( nbdHist->FindBin(i), nbdDist->Eval(i) );
    }
 
  return nbdHist;

}
