#!/bin/bash

#Runs CondorSubmit.xml and converts picoDsts from catalog query to femtoDsts

config_file=../../femto.config
############################################################################
#RUNNING WITH CONFIG FILE
############################################################################

source $config_file

MultDir=$PARENT_DIR/multiplicity_histograms

if [ ! -d $MultDir ]; then
    echo "Creating the mult directory $ProfileDir"
    mkdir $MultDir
fi

outDir=$MultDir/out/
logDir=$MultDir/log/
schedDir=$MultDir/sched/
fDstDir=$FDST_OUT 
workDir=$ANA_DIR
inputParameter=$workDir/input_parameters/parameter_n0p5_0_norm.list

filelist=femto_3GeV.list
#Check to make sure that the inputFileList exists 
if [ ! -e $filelist ]; then
    readlink -f ${fDstDir}/*.root > femto_3GeV.list

fi

filelist=`readlink -f femto_3GeV.list`
echo $filelist

packageName="multDists"

#Check to Make sure the output Directory exists
if [ ! -d $outDir ]; then
    echo "Creating the out directory $outDir"
    mkdir $outDir
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

#echo $packageName
star-submit-template -template ../condor_xml/MultSubmitter.xml -entities logDir=$logDir,filelist=$filelist,outDir=$outDir,schedDir=$schedDir,inputParameter=$inputParameter,packageName=$packageName,workDir=$workDir
