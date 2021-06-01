#!/bin/bash

#Runs CondorSubmit.xml and converts picoDsts from catalog query to femtoDsts
#If running from config file, just run ./RunCondorSubmit.sh. 

#If specifying the directories, fill the parameters as follows: 

outDir="fDsts output directory ... $1"
logDir="fDsts log directory ... $2"
schedDir="fDsts schedular directory ... $3"
goodRunList="fDsts good run file list ... $4"
goodRunTxt="fDsts good run numbers ... $5"

if [ "$#" -ne 0 ] && [ "$#" -ne 5 ]; then
    echo "Bad number of arguments..." 
    echo "Either run from config (no arguments)"
    echo "or specify all five parameters"
    exit
fi 

config_file=../femto.config
############################################################################
#RUNNING WITH CONFIG FILE
############################################################################
if [ "$#" -eq 0 ]; then

    source $config_file
    fDstDir=$PARENT_DIR/fDst

    if [ ! -d $fDstDir ]; then
	echo "Creating the log directory $fDstDir"
	mkdir $fDstDir
    fi

    outDir=$fDstDir/out
    logDir=$fDstDir/log
    schedDir=$fDstDir/sched
    goodRunList=$GOOD_RUNS_LIST
    goodRunTxt=$GOOD_RUNS_TXT

fi

#RUNNING WITH SPECIFIED DIRECTORIES
###########################################################################
#Check to make sure there are the right number of arguments
###########################################################################
if [ "$#" -eq 5 ]; then

    outDir=$1
    logDir=$2
    schedDir=$3
    goodRunList=$4 
    goodRunTxt=$5 # probably should be filelist/GoodRuns.txt
fi    

#Check to Make sure the output Directory exists
if [ ! -d $outDir ]; then
    echo "Creating the log directory $outDir"
    mkdir $outDir
fi

#Check to make sure that the inputFileList exists 
if [ ! -e $goodRunList ]; then
    echo "ERROR: Input file list $goodRunList does not exist!"
    exit 1
fi

#Check to make sure that the inputFileList exists 
if [ ! -e $goodRunTxt ]; then
    echo "ERROR: Input file list $goodRunTxt does not exist!"
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

echo "      Out directory: " $outDir
echo "      Log directory: " $logDir
echo "   Sched. directory: " $schedDir
echo "       Good RunList: " $goodRunList
echo "    Good RunNumbers: " $goodRunTxt

#export the directories to the config file
echo "" >> $config_file
echo "# Generated from RunCondorSubmit.sh" >> $config_file
echo "FDST_OUT="$outDir >> $config_file


#star-submit-template -template CondorSubmit.xml -entities logDir=$logDir,outDir=$outDir,schedDir=$schedDir,goodRunList=$goodRunList,goodRunTxt=$goodRunTxt
