#!/bin/bash

outDir=/star/data01/pwg/sheppel/femtoAnalysis/out_sys_5_18/
logDir=/star/data01/pwg/sheppel/femtoAnalysis/log/
schedDir=/star/data01/pwg/sheppel/femtoAnalysis/sched/
filelist=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/input_parameters/femto_3GeV.list

configsArray=("n0p2_0" "n0p3_0" "n0p4_0" "n0p5_0" "n0p5_0_pt1" "n0p5_0_pt2" "n0p5_0_pt3" "n0p5_0_pt4")

for configs in ${configsArray}; do
    ./sys_submitter.sh $outDir $logDir $schedDir $filelist ${configs} 
    wait
done

