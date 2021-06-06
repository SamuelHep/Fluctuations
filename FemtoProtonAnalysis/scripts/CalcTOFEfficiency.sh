#!/bin/bash

###Sets up the tofmatch efficiency file####

config_file=../../femto.config
source $config_file

tofmatch_dir=$TOFMATCH_DIR/out/

if [ ! -d $tofmatch_dir ]; then
    echo "Directory doesnt exist, did you run ./run_tof_submitter.sh?"
fi

cp MadAdder.py $tofmatch_dir/.
current_dir=`pwd`
cd $tofmatch_dir
python MadAdder.py tofMatch_ tofmatch_efficiency.root
wait
cp tofmatch_efficiency.root $ANA_DIR/rootfiles/eff/.
cd $current_dir





