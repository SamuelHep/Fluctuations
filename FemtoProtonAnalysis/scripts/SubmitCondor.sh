#!/bin/bash

#Runs CondorSubmit.xml and converts picoDsts from catalog query to femtoDsts

outDir=$1
logDir=$2
schedDir=$3
filelist=$4 
inputParameter=$5
packageName=$6

#Check to make sure there are the right number of arguments
if [ "$#" -ne 6 ]; then
    echo "ERROR: Incorrect number of arguments! This script needs four: outputDir logDir schedDir filelist inputParameter packageName"
    exit 1
fi

#Check to Make sure the output Directory exists
if [ ! -d $outDir ]; then
    echo "ERROR: Output Directory $outDir does not exist!"
    exit 1
fi

#Check to make sure that the inputFileList exists 
if [ ! -e $filelist ]; then
    echo "ERROR: Input file list $filelist does not exist!"
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

#Check that the scheduler directory exists, if not make it
if [ ! -e $inputParameter ]; then
    echo "ERROR: inputParameter does not exist!"
    exit 1
fi

star-submit-template -template ../condor_xml/CondorSubmit.xml -entities logDir=$logDir,filelist=$filelist,outDir=$outDir,schedDir=$schedDir,inputParameter=$inputParameter,packageName=$inputPackageName
