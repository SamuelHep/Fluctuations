#!/bin/bash

#If running from the femto.config file, just run ./submit_full_analysis.sh

outDir="Out directory of moment profiles ... $1"
logDir="log directory ... $2"
schedDir="scheduler directory ... $3"
fDstDir="directory with fDsts rootfiles ... $4"
tpc_eff_file="TPC efficiency file ... $5"
tof_eff_file="TOF efficiency file ... $6"

if [ "$#" -ne 0 ] && [ "$#" -ne 7 ]; then
    echo "Bad number of arguments..." 
    echo "Either run from config (no arguments)"
    echo "or specify all five parameters"
    exit
fi 

config_file=../../femto.config
############################################################################
#RUNNING WITH CONFIG FILE
############################################################################
if [ "$#" -eq 0 ]; then

    source $config_file

    ProfileDir=$PARENT_DIR/moment_profiles

    if [ ! -d $ProfileDir ]; then
	echo "Creating the log directory $ProfileDir"
	mkdir $ProfileDir
    fi

    outDir=$ProfileDir/out/
    logDir=$ProfileDir/log/
    schedDir=$ProfileDir/sched/
    fDstDir=$FDST_OUT 
    tpc_eff_file=$TPC_EFF_FILE
    tof_eff_file=$TOF_EFF_FILE

fi

#RUNNING WITH SPECIFIED DIRECTORIES
###########################################################################
#Check to make sure there are the right number of arguments
###########################################################################
if [ "$#" -eq 7 ]; then

    source $config_file

    outDir=$1
    logDir=$2
    schedDir=$3
    fDstDir=$4
    tpc_eff_file=$5
    tof_eff_file=$6

fi

filelist=femto_3GeV.list
#Check to make sure that the inputFileList exists 
if [ ! -e $filelist ]; then
    readlink -f ${fDstDir}/*.root > femto_3GeV.list

fi

filelist=`readlink -f femto_3GeV.list`
echo $filelist

configsArray=("n0p2_0" "n0p3_0" "n0p4_0" "n0p5_0" "n0p5_0_pt1" "n0p5_0_pt2" "n0p5_0_pt3" "n0p5_0_pt4")

for configs in ${configsArray}; do
    ./sys_submitter.sh $outDir $logDir $schedDir $filelist $tpc_eff_file $tof_eff_file ${configs} 
    wait
done
