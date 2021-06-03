#!/bin/bash

###Takes the reduced embedding files and makes a proton efficency file###

./hadd_embedding.sh

config_file=../../femto.config
source $config_file

embedding_file=$EMBED_TREE

setup 64bits

root -l -q -b ../macros/RunProcessEmbedding.C\(\"$embedding_file\",\"$ANA_DIR\"\)
wait

hadd -f ../rootfiles/eff/tpc_efficiency.root ../rootfiles/eff/protonEmbedding_nhf10.root ../rootfiles/eff/protonEmbedding_nhf12.root ../rootfiles/eff/protonEmbedding_nhf15.root


