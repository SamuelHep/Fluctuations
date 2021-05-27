#!/bin/bash

outDir=/star/data01/pwg/sheppel/femtoAnalysis/out_additional_tests/
logDir=/star/data01/pwg/sheppel/femtoAnalysis/log/
schedDir=/star/data01/pwg/sheppel/femtoAnalysis/sched/
filelist=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/femto_3GeV.list
inputParameter=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p9_n0p5_TPC.list
packageName=n0p9_n0p5_TPC

./JobSubmitter.sh $outDir $logDir $schedDir $filelist $inputParameter $packageName
