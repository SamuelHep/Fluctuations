#!/bin/bash

outDir=/star/data01/pwg/sheppel/femtoAnalysis/out_-0p2_0p1_norm/
logDir=/star/data01/pwg/sheppel/femtoAnalysis/log/
schedDir=/star/data01/pwg/sheppel/femtoAnalysis/sched/
filelist=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/femto_3GeV.list
inputParameter=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p2_0p1_norm.list
packageName=norm_n0p2_0p1

./JobSubmitter.sh $outDir $logDir $schedDir $filelist $inputParameter $packageName
