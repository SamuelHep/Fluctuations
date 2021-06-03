#!/bin/bash

setup 64bits

./hadd_profiles.sh

config_file=../../femto.config
source $config_file

$profile_dir_local=$ANA_DIR/rootfiles/profiles
$cumulant_dir_local=$ANA_DIR/rootfiles/cumulants
$pileup_corr=$PU_NORM_FILE
$pileup_corr_high=$PU_HIGH_FILE
$pileup_corr_low=$PU_LOW_FILE

root -l -q -b ../macros/RunPileUpCorrections.C\(\"$profile_dir_local\", \"$cumulant_dir_local\", \"$pileup_corr\", \"$pileup_corr_low\", \"$pileup_corr_high\"\)
wait

root -l -q -b ../macros/CalcSystematic.C\(\"$cumulant_dir_local\"\)
wait

root -l -q -b ../macros/MakeTopFive.C
wait

root -l -q -b ../macros/DrawCumulant.C


