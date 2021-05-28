#!/bin/bash

outDir=$1
logDir=$2
schedDir=$3
filelist=$4
name=$5
workDir=`pwd`/..
inputParameter0=$workDir/input_parameters/parameter_${name}_norm.list
inputParameter1=$workDir/input_parameters/parameter_${name}_SYS1.list
inputParameter2=$workDir/input_parameters/parameter_${name}_SYS2.list
inputParameter3=$workDir/input_parameters/parameter_${name}_SYS3.list
inputParameter4=$workDir/input_parameters/parameter_${name}_SYS4.list
inputParameter5=$workDir/input_parameters/parameter_${name}_SYS5.list
inputParameter6=$workDir/input_parameters/parameter_${name}_SYS6.list
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
echo " "
echo " package name: " $packageName
echo " "

sleep 4s

echo "Submitting..."

./SysSubmitter.sh $outDir $logDir $schedDir $filelist $packageName $inputParameter0 $inputParameter1 $inputParameter2 $inputParameter3 $inputParameter4 $inputParameter5 $inputParameter6
