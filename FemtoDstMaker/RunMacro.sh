#!/bin/bash
#simple test script for generating femtoDsts
#filelist=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/3GeV_newProd_Fluct_GoodList.list
#runfile=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/good_3GeV.txt

filelist=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/test.list
runfile=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/test.txt


nEvents=99999999  # number of events 
outdir=/star/data01/pwg/sheppel/femtoDsts_2_15/      # out directory
ID=Fxt_3GeV         # output file label

root4star -q -b StFemtoMacro.C\($nEvents,\"$filelist\",\"$outdir/\",\"$ID\",\"$runfile\",\"checkEmbedding_pip.root\"\) 


