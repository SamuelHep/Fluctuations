#!/bin/bash

p0w5=../rootfiles/profiles/smallSample_profiles.root 
o_p0w5_tm=../rootfiles/cumulants/smallSample_out.root

toyHists_nocut=../rootfiles/pileupCor/CorPlots_nocut.root

root -l -q -b ../macros/RunPileUpCorrection_bs.C\(\"$p0w5\",\"$o_p0w5_tm\",\"$toyHists_nocut\",0\)

