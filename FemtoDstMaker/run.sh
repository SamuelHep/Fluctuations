#!/bin/bash

outDir=/star/data01/pwg/sheppel/femtoDsts_2_15/      # out directory
logDir=/star/data01/pwg/sheppel/femtoDsts_11_9_log/
schedDir=/star/data01/pwg/sheppel/femtoDsts_11_9_sched/
goodRunList=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/good_3GeV.list
goodRunTxt=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/good_3GeV.txt

./RunCondorSubmit.sh $outDir $logDir $schedDir $goodRunList $goodRunTxt
