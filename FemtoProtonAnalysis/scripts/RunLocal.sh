#!/bin/bash

FILELIST=test.list
parameterList=../input_parameters/parameter_n0p5_0_norm.list
infilename=/star/data01/pwg/sheppel/codeTest/fDst/out/594C9CE754246CDBDA770F84C8455D42_351.fDst.root
tpc_eff=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/eff/tpc_efficiency.root
tof_eff=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/eff/tof_efficiency_fineBinning.root

root -b -l -q ../macros/RunEvents.C\(\"$FILELIST\",\"$parameterList\",\"$infilename\"\,\"$tpc_eff\",\"$tof_eff\"\)
