#!/bin/bash

p0w5=rootfiles/profiles_0w5.root
o_p0w5_tm=rootfiles/cumulants_0w5_toyModel_gmOnly.root

toyHists_nocut=rootfiles/CorPlots_nocut_gmOnly.root
root -l -q -b ../macros/RunPileUpCorrection_bs.C\(\"$p0w5\",\"$o_p0w5_tm\",\"$toyHists_nocut\",0\)

