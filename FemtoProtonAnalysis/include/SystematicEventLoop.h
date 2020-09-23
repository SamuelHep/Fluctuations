#ifndef SYS_EVENTLOOP_H
#define SYS_EVENTLOOP_H

#include "TObject.h"
#include "InputParameterList.h"
#include <vector>
#include "TString.h"

class TChain;
class CumulantProfileContainer;
class ProtonEfficiency;
class TRandom3;


void GetProfilesAndCalculateSystematics(TString infilename,TString outfilename, int nBootstraps, TString sysLabel, double start, double sysDelta);

int SystematicEventLoopPrintToFile(
				   TChain * tc,
				   long int nentries,
				   InputParameterList & pl,
				   ProtonEfficiency * eff,
				   TRandom3 * rand,
				   int nBootstraps,
				   TString sysLabel,
				   TString outfileName
				   );


int SysEventLoop(
		 TChain * tc,
		 long int nentries,
		 InputParameterList & pl,
		 std::vector<CumulantProfileContainer*> cpc_vec,
		 int nBootstraps,
		 TString sysLabel,
		 ProtonEfficiency * eff,
		 TRandom3 * rand
		 );


#endif
