#!/bin/bash

outDir=/star/data01/pwg/sheppel/femtoAnalysis/out/
logDir=/star/data01/pwg/sheppel/femtoAnalysis/log/
schedDir=/star/data01/pwg/sheppel/femtoAnalysis/sched/
filelist=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/femto_3GeV.list
inputParameter0=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p6_n0p4_norm.list
inputParameter1=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p6_n0p4_SYS1.list
inputParameter2=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p6_n0p4_SYS2.list
inputParameter3=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p6_n0p4_SYS3.list
inputParameter4=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p6_n0p4_SYS4.list
inputParameter5=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/parameter_n0p6_n0p4_SYS5.list
packageName=norm_n0p6_n0p4

./SysSubmitter.sh $outDir $logDir $schedDir $filelist $packageName $inputParameter0 $inputParameter1 $inputParameter2 $inputParameter3 $inputParameter4 $inputParameter5
