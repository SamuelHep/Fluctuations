#!/bin/bash
#Submit for an analysis window and systematics
 
outDir=$1
logDir=$2
schedDir=$3
filelist=$4
tpc_eff_file=$5
tof_eff_file=$6
name=$7
workDir=`pwd`/..
inputParameter0=$workDir/input_parameters/parameter_${name}_norm.list
inputParameter1=$workDir/input_parameters/parameter_${name}_SYS1.list
inputParameter2=$workDir/input_parameters/parameter_${name}_SYS2.list
inputParameter3=$workDir/input_parameters/parameter_${name}_SYS3.list
inputParameter4=$workDir/input_parameters/parameter_${name}_SYS4.list
inputParameter5=$workDir/input_parameters/parameter_${name}_SYS5.list
inputParameter6=$workDir/input_parameters/parameter_${name}_SYS6.list
inputParameter7=$workDir/input_parameters/parameter_${name}_SYS7.list
inputParameter8=$workDir/input_parameters/parameter_${name}_SYS8.list
packageName=${name}

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Out   directory: " $outDir   
echo "Log   directory: " $logDir   
echo "Sched directory: " $schedDir   
echo " "
echo "FileList = " $filelist
echo " "
echo "Input files:" 
echo "     " $inputParameter0
echo "     " $inputParameter1
echo "     " $inputParameter2
echo "     " $inputParameter3
echo "     " $inputParameter4
echo "     " $inputParameter5
echo "     " $inputParameter6
echo "     " $inputParameter7
echo "     " $inputParameter8
echo " "
echo " package name: " $packageName
echo " "

sleep 4s

#Check to make sure there are the right number of arguments
if [ "$#" -ne 7 ]; then
    echo "ERROR: Incorrect number of arguments! This script needs five: outputDir logDir schedDir filelist packageName"
    exit 1
fi

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

#Check to make sure that the tpc efficiency exists 
if [ ! -e $tpc_eff_file ]; then
    echo "ERROR: Input file list $tpc_eff_file does not exist!"
    exit 1
fi

#Check to make sure that the tof efficiency exists 
if [ ! -e $tof_eff_file ]; then
    echo "ERROR: Input file list $tof_eff_file does not exist!"
    exit 1
fi

star-submit-template -template ../condor_xml/SysJobSubmitter.xml -entities logDir=$logDir,filelist=$filelist,outDir=$outDir,schedDir=$schedDir,packageName=$packageName,workDir=$workDir,tpc_eff_file=$tpc_eff_file,tof_eff_file=$tof_eff_file,inputParameter0=$inputParameter0,inputParameter1=$inputParameter1,inputParameter2=$inputParameter2,inputParameter3=$inputParameter3,inputParameter4=$inputParameter4,inputParameter5=$inputParameter5,inputParameter6=$inputParameter6,inputParameter7=$inputParameter7,inputParameter8=$inputParameter8

