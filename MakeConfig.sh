#!/bin/bash
#This script generates a config file that is used by both the FemtoDstMaker and FemtoProtonAnalysis.
#It will also setup the enviroment and compile all the code

#### EDIT BEFORE RUNNING ###########################################################

#Make this a location with ~200 Gb ( like /star/data01/pwg/your_name/femtoAnalysis/ )
output_parent_directory=/star/data01/pwg/sheppel/test_femtoAnalysis

# External files needed for running (probably don't need to change):
tpc_efficiency_file=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/eff/tpc_efficiency.root
tof_efficiency_file=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/eff/tof_efficiency_fineBinning.root
pileup_correction_file=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/pileupCor/test_data_best.root
pileup_correction_file_syslow=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/pileupCor/test_data_best_minus.root
pileup_correction_file_syshigh=/star/u/sheppel/femtoRepo/FemtoProtonAnalysis/rootfiles/pileupCor/test_data_best_plus.root
good_run_file_list=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/3GeV_newProd_Fluct_GoodList.list
good_runs_txt=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/good_3GeV.txt

# If you can't find the efficiency files or the pileup_correction_files, they can be located here:
# Efficiency files -> https://drupal.star.bnl.gov/STAR/system/files/EfficiencyCorr.zip
# PileUp Correction files -> https://drupal.star.bnl.gov/STAR/system/files/pileupCorr.zip
# GoodRun files -> https://drupal.star.bnl.gov/STAR/system/files/good_run_lists.zip
####################################################################################

nocompile=false 
makeclean=false
while [ -n "$1" ]; do 
    case "$1" in
	-nc) nocompile=true ;;
	*) echo "Option $1 not recognized"
    esac
    shift 
done
 
### Lets make sure we have everything
ext_files=($tpc_efficiency_file 
$tof_efficiency_file
$pileup_correction_file
$pileup_correction_file_syslow
$pileup_correction_file_syshigh
$good_run_file_list
$good_runs_txt )
all_set=true

for file in ${ext_files[@]}; do
    if [ ! -f $file ]; then
	echo "File Not found! Can't find: " $file
	all_set=false
    fi
done

if [ $all_set = false ]; then
    echo "Please locate missing file/files and run again... EXITING"
    exit
fi
 
### Make sure the output directory exists or it can be created ###

if [ ! -d $output_parent_directory ]; then
    echo "Creating output_parent_directory = $output_parent_directory"
    mkdir $output_parent_directory
fi

if [ ! -d $output_parent_directory ]; then
    echo "Creating directory failed... EXITING"
    exit
fi

top_dir=`pwd`
maker_dir=${top_dir}/FemtoDstMaker/
analysis_dir=${top_dir}/FemtoProtonAnalysis/

# Generate config file
config_file=femto.config

if [ -f $config_file ]; then
    rm $config_file
fi

touch $config_file
echo "# Generated from MakeConfig.sh" >> $config_file
echo "TOP_DIR="$top_dir >> $config_file
echo "MAKER_DIR="$maker_dir >> $config_file
echo "ANA_DIR="$analysis_dir >> $config_file
echo "PARENT_DIR="$output_parent_directory >> $config_file
echo "TPC_EFF_FILE="$tpc_efficiency_file >> $config_file
echo "TOF_EFF_FILE="$tof_efficiency_file >> $config_file
echo "PU_NORM_FILE="$pileup_correction_file >> $config_file
echo "PU_HIGH_FILE="$pileup_correction_file_syshigh >> $config_file
echo "PU_LOW_FILE="$pileup_correction_file_syslow >> $config_file
echo "GOOD_RUNS_LIST="$good_run_file_list >> $config_file
echo "GOOD_RUNS_TXT="$good_runs_txt >> $config_file

# Change to 64bits #
setup 64bits

if [ $nocompile = true ]; then
    echo "Not compiling and exiting early."
    exit
fi

# Compile the FemtoMaker #
cd $maker_dir
cons

# Compile the Analysis code #
cd $analysis_dir
make

#back home
cd $top_dir
