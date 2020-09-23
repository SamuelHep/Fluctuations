#ifndef EVENTLOOP_H
#define EVENTLOOP_H

#include "TObject.h"
#include "InputParameterList.h"
#include <vector>

class TChain;
class CumulantProfileContainer;
class ProtonEfficiency;
class TRandom3;

int EventLoop(
	      TChain * tc,
	      long int nentries,
	      InputParameterList & pl,
	      CumulantProfileContainer * cpc,
	      ProtonEfficiency * eff
	      );


int EventLoopSystematic(
			TChain * tc,
			long int nentries,
			InputParameterList & pl,
			vector<CumulantProfileContainer*> &cpc,
			ProtonEfficiency * eff
			);


int EventLoopBootstrap(
		       TChain * tc,
		       long int nentries,
		       InputParameterList & pl,
		       CumulantProfileContainer* cpc,
		       std::vector<CumulantProfileContainer*> cpc_vec,
		       ProtonEfficiency * eff,
		       TRandom3 * rand
		       );

int EventLoopSimPoisson(
			long int nentries_perCentBin,
			CumulantProfileContainer * cpc,
			TRandom3 * rand,
			int minCent,
			int maxCent,
			double eff
			);


#endif
