# Pile Up Correction File #

## Requires the FxtMult.root and an input distribution ##

FxtMult.root should be generated from the FemtoDsts by running 
./run_mult_submitter.sh in FemtoProtonAnalysis and ./hadd_mult.sh.
Copy FxtMult.root into this directory.

The input distribution can be any probability distribution like a uniform
distribution. /star/u/glauber_pu1per_eff.root will work.

This macro performs unfolding to extract true refmult3 distribution 
from sinle-collision events, based on experimental measured distribution 
including pilupe events.

  usage:
  root -l test_data_best.C++

  for error estimation:
  root -l test_data_best.C++'(1,iBS)'
  , where iBS is integer number from 0-99.

One can repeat the unfolding by varying pileup probability, alpha, to 
determine the best parameter which can describe the data.

