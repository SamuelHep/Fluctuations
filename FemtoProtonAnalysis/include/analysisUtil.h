#ifndef ANALYSIS_UTIL_H
#define ANALYSIS_UTIL_H

#include <vector>
#include "TChain.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "CumulantProfileContainer.h"
#include "InputParameterList.h"

float getProtonDedxMeanShift(float p);

TChain * GetTChainFromList(TString listname, TString treename);

int ComputeBootstrapErrors(CumulantProfileContainer * primaryCpc, std::vector<CumulantProfileContainer*> cpc_vec,
			   std::vector<TGraphErrors*> & compare_graphs);

int BootstrapGraph(TGraphErrors * gr, std::vector<TGraphErrors*> gr_vec);

int BootstrapGraph(TGraphErrors * gr, std::vector<TGraphErrors*> gr_vec, TGraphErrors * gr_err, std::vector<TGraphErrors*> & compare_graphs);

std::vector<TGraphErrors*> MakeSysComparisonPlots(
					     CumulantProfileContainer * primary,
					     std::vector<CumulantProfileContainer*> sys,
					     TString sysLabel,
					     double start,
					     double sysDelta
						  );

InputParameterList ReadInputFile(TString inputfile);

void RunPileUpCorr(TString infilename, TString outfilename, TString histfilename,Bool_t urqmdHists);

int GetRunIndex( std::vector<int> &runvec, int runNum);
std::pair<bool,bool> PileUpBadEvent(int fxt,double nmip,double pipdu);
std::pair<bool,bool> PileUpBadEvent(int gmult, int fxt,double nmip,double pipdu);
std::pair<bool,bool> PileUpBadEventVariable(int fxt,double nmip,double pipdu,int index);

int CentBin(int mult);
int CentBin3(int mult);
int CentBinTofMatch(int mult);


#endif
