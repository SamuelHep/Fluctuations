#!/bin/bash

#Runs CondorSubmit.xml and converts picoDsts from catalog query to femtoDsts

outDir=$1
logDir=$2
schedDir=$3
goodRunList=$4 # probably should be filelist/GoodRuns.txt

#Check to make sure there are the right number of arguments
if [ "$#" -ne 4 ]; then
    echo "ERROR: Incorrect number of arguments! This script needs four: outputDir logDir schedDir goodRunsList"
    exit 1
fi

#Check to Make sure the output Directory exists
if [ ! -d $outDir ]; then
    echo "ERROR: Output Directory $outDir does not exist!"
    exit 1
fi

#Check to make sure that the inputFileList exists 
if [ ! -e $goodRunList ]; then
    echo "ERROR: Input file list $goodRunList does not exist!"
    exit 1
fi

#Check that the log directory exists, if not make it
if [ ! -d $logDir ]; then
    echo "Creating the log directory $logDir"
    mkdir $logDir
fi

#Check that the scheduler directory exists, if not make it
if [ ! -d $schedDir ]; then
    echo "Creating the scheduler directory $schedDir"
    mkdir $schedDir
fi

star-submit-template -template CondorSubmit.xml -entities logDir=$logDir,goodRunList=$goodRunList,outDir=$outDir,schedDir=$schedDir
