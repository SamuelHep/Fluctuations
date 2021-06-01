#!/bin/bash

outDir=/star/data01/pwg/sheppel/femtoAnalysis/tof_eff/out_eta/
logDir=/star/data01/pwg/sheppel/femtoAnalysis/tof_eff/log/
schedDir=/star/data01/pwg/sheppel/femtoAnalysis/tof_eff/sched/
filelist=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/femto_3GeV.list
inputParameter=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p5_0_norm.list
packageName=tof_eff

./TOFSubmitter.sh $outDir $logDir $schedDir $filelist $inputParameter $packageName
