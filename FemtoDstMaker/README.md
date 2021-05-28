
# FemtoDstMaker
## generates a reduced file format for the proton moments analysis

The code includes a few makers:
1. StFemtoDstMaker     
*takes picoDst and generates reduced rootfile format (fDst), a tree of StFemtoEvents*
2. StFemtoEvent        
*reduced event class, Event info and array of StFemtoTracks*
4. StFemtoTrack
*reduced track class, mainly momentum and dE/dx nSigmas*
6. EpdPileUpRejection  
*just used to calculate nMipEPD. Not used in analysis*
8. StEpdUtil
*just used to calculate nMipEPD. Not used in analysis*

## How to compile and run the code:
  
First, the code needs to be compiled. Run the following lines:

`setup 64bits`

`cons`

Second, edit the run.sh and change the file paths to convient out, log and scheduler directories (outDir, logDir and schedDir):

> outDir=/star/data01/pwg/your_directory/fDsts/  
> logDir=/star/data01/pwg/your_directory/fDsts_log/
> schedDir=/star/data01/pwg/your_directory/fDsts_sched/
> goodRunList=/star/u/sheppel/femtoRepo/FemtoDstMaker/filelist/good_3GeV.txt

You can leave the goodRunList as the current path or obtain the good run list from https://drupal.star.bnl.gov/STAR/system/files/good_3GeV.txt

Then, it you just run ./run.sh and wait for the jobs to finish

`./run.sh`

After a while, 385 rootfiles should be produced.
