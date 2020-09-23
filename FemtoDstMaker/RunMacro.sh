#!/bin/bash
#simple test script for generating femtoDsts

filelist=filelist/test.list
runfile=filelist/test.fxt

nEvents=999999  # number of events 
outdir=./       # out directory
ID=test         # output file label

root4star -q -b StFemtoMacro.C\($nEvents,\"$filelist\",\"$outdir/\",\"$ID\",\"$runfile\"\) 


