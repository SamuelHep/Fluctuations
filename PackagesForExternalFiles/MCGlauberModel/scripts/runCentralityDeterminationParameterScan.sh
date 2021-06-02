#!/bin/bash

#This script is used to find the centrality binning of data
#It needs a Glauber Simulation to have been completed first.

#USER DEFINED VARIABLES
source ../../../femto.config

dataFileName=$FXTMULT3_FILE
dataHistoName=fxtmult3       #Name of the histogram in the file from above
glauberFileName=../data/Glauber_197_197_27.3mb_WoodsSaxon.root #Glauber file generated from runGlauberSimulation.sh
outputFileName=../data/GlauberFit.root

normStartBinCenter=20                                      #The bin center value to begin the chi^2 matching/optimization routine
normStopBinCenter=80                                      #The bin center value to end the chi^2 matching/optimization routine. Use -1 to use the last bin of the data histogram. NOTE: This value must be larger than normStartBinCenter to make any sense.
numberOfPointsToEvaluate=30                               #The number of n,k pairs to analyze. The default is 500               
#yourEmail=zacharysweger@gmail.com

npp=0.6
nppLow=0.4
nppHigh=0.8

k=0.1
kLow=0.0
kHigh=10.0

collisonhardness=0.06
collisionhardnessLow=0.0
collisionhardnessHigh=1.0

freeMu=true
freeK=true
freeX=false

doCentralityCuts=true
doTheFullAnalysis=false

outputFileName=OFILE.root
root -l -b -q ../macros/RunCentralityDetermination.C\(\"$dataFileName\",\"$dataHistoName\",\"$glauberFileName\",\"$outputFileName\",$normStartBinCenter,$normStopBinCenter,$numberOfPointsToEvaluate,$npp,$k,$collisonhardness,$nppLow,$nppHigh,$kLow,$kHigh,$collisionhardnessLow,$collisionhardnessHigh,$doTheFullAnalysis,$doCentralityCuts,$freeMu,$freeK,$freeX\)


