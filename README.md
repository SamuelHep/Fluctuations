# Fluctuations
STAR framework code to make .femto.root files and preform proton fluctuation analysis 

## FemtoDstMaker 
Creates a reduced file format for the proton analysis.

## FemtoProtonAnaylsis
Proton Fluctuations analysis for FXT 3 GeV

## BEFORE RUNNING ANYTHING
Edit the `MakeConfig.sh` file. 
An output directory must be specified.
The external files may need to be changed. If so, the location of the files are stated in the script.

## How to run.... ##

First, run in your terminal:

> setup 64bits 

You must do this every session to run anything in FemtoProtonAnalysis.
If you compile and do not do this, it will fail and you need to go into `FemtoProtonAnaylsis` and `make clean`

Next, edit the paths in `MakeConfig.sh`. Change the output_parent_directory to somewhere with space ~200 Gb.
You can leave the external files as is.

Run `./MakeConfig.sh` to format the *femto.config* file.

Go into *FemtoDstMaker* and run:

> ./RunCondorSubmit.sh 

wait a few hours for the reduced file format to be generated.

Go into *FemtoProtonAnalysis/scripts* and run:

> ./submit_full_analysis.sh

wait a few hours for the moment profiles to be generated.

Next, run the code that adds the profiles to a single root file, runs the corrections (pileup and CBWC) and plots the results:

> ./run_corr_and_plot.sh 

A plot with the rapidity scan and pT scan is generated in `FemtoProtonAnaylsis/img/`

For rootfiles with the results, look at cumulant*.root files in `FemtoProtonAnaylsis/rootfiles/cumulants/`

## External Packages ##

To generate the centrality definition, efficiency curves and pileup correction file, you will need to use the external packages.
