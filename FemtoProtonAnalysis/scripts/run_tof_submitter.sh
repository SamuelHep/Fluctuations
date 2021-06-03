#!/bin/bash


config_file=../../femto.config
############################################################################
#RUNNING WITH CONFIG FILE
############################################################################
if [ "$#" -eq 0 ]; then

    source $config_file

    tofMatchDir=$PARENT_DIR/tofmatch

    if [ ! -d $tofMatchDir ]; then
	echo "Creating the log directory $tofMatchDir"
	mkdir $tofMatchDir
    fi

    outDir=$tofMatchDir/out
    logDir=$tofMatchDir/log
    schedDir=$tofMatchDir/sched
    fDstDir=$FDST_OUT 

fi

filelist=femto_3GeV.list
#Check to make sure that the inputFileList exists 
if [ ! -e $filelist ]; then
    readlink -f ${fDstDir}/*.root > femto_3GeV.list

fi

inputParameter=$ANA_DIR/input_parameters/parameter_n0p5_0_norm.list
packageName=tof_eff


#Check to Make sure the output Directory exists
if [ ! -d $outDir ]; then
    echo "Creating the log directory $outDir"
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

workDir=$ANA_DIR

star-submit-template -template ../condor_xml/TOFEfficiencyMaker.xml -entities logDir=$logDir,filelist=$filelist,outDir=$outDir,schedDir=$schedDir,inputParameter=$inputParameter,packageName=$packageName,workDir=$workDir
