#ifndef MOMENTEVENTLOOP_H
#define MOMENTEVENTLOOP_H

#include "TObject.h"
#include "InputParameterList.h"
#include <vector>
#include "TString.h"

class TChain;
class CumulantProfileContainer;
class ProtonEfficiency;
class TRandom3;

int MomentEventLoopLocal(
			 TChain * tc,
			 long int nentries,
			 InputParameterList & pl,
			 ProtonEfficiency * eff,
			 TRandom3 * rand,
			 int nBootstraps,
			 TString outfileName
			 );

int MomentEventLoop(
		    TChain * tc,
		    long int nentries,
		    InputParameterList & pl,
		    CumulantProfileContainer* cpc,
		    std::vector<CumulantProfileContainer*> cpc_vec,
		    ProtonEfficiency * eff,
		    TRandom3 * rand
		    );


int MomentEventLoopPrintToFile(
			       TChain * tc,
			       long int nentries,
			       InputParameterList & pl,
			       ProtonEfficiency * eff,
			       TRandom3 * rand,
			       int nBootstraps,
			       TString outfileName
			       );


int PileUpEventLoop(
		    TChain * tc,
		    long int nentries,
		    InputParameterList & pl,
		    TString outfilename
		    );


#endif
