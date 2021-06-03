#!/bin/bash

config_file=../../femto.config
source $config_file
currentDir=`pwd`
$profile_dir=$PROFILE_DIR_OUT
$profile_dir_local=$ANA_DIR/rootfiles/profiles

cp AddAll.sh $profile_dir
cp MadAdder.y $profile_dir

cd $profile_dir
./AddAll.sh
wait

cp profiles_* $profile_dir_local
cd $currentDir
