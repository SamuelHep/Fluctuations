#!/bin/bash

config_file=../../femto.config
source $config_file
currentDir=`pwd`

MultDir=$PARENT_DIR/multiplicity_histograms/out/

cp MadAdder.py $MultDir

cd $MultDir

python MadAdder.py FxtMult FxtMult.root
wait

mult_hist=`readlink -f FxtMult.root`
cd $currentDir

echo "" >> $config_file
echo "#Mult hist generated" >> $config_file
echo "FXTMULT3_FILE="$mult_hist >> $config_file

