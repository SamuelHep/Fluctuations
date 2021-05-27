#ifndef QALOOP_H
#define QALOOP_H

#include "TObject.h"
#include "InputParameterList.h"
#include <vector>
#include <utility>

class TChain;

int SimplePlots(
		TChain * tc,
		long int nentries,
		InputParameterList & pl,
		TString outfilename
		);

int QALoop(
	   TChain * tc,
	   long int nentries,
	   InputParameterList & pl,
	   TFile * outfile
	   );

int QAplots(
	    TChain * tc,
	    long int nentries,
	    InputParameterList & pl,
	    TString outfilename
	    );


#endif
