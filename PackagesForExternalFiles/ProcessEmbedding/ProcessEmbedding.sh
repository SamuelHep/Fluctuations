#!/bin/bash

config_file=../../femto.config
source $config_file

PD=$EMBED_DIR

outDir=$PD/out/
logDir=$PD/log/
schedDir=$PD/sched/
workDir=`pwd`

if [ ! -d $PD ]; then
    echo "Creating the parent directory $PD"
    mkdir $PD
fi

#Check to Make sure the output Directory exists
if [ ! -d $outDir ]; then
    echo "Creating the log directory $outDir"
    mkdir $outDir
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

star-submit-template -template muEmbed.xml -entities outDir=$outDir,logDir=$logDir,schedDir=$schedDir,workDir=$workDir 
