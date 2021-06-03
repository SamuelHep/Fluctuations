#!/bin/bash

config_file=../../femto.config
source $config_file
currentDir=`pwd`
embed_dir=$EMBED_DIR/out/

cp MadAdder.py $embed_dir

cd $embed_dir

python MadAdder.py muEmbed_ muEmbed_proton.root
wait

embed_tree=`readlink -f muEmbed_proton.root`
cd $currentDir

echo "" >> $config_file
echo "#Reduced proton embedding tree generated" >> $config_file
echo "EMBED_TREE="$embed_tree >> $config_file

