
# FemtoProtonAnalysis #

## Calculates the proton cumulants from fDst reduced file format ##

General code flow:

1. Generate a file list of fDsts
2. Run through all events and generate profiles for each averaged quantity i.e. <m_1_1>, ..., <m_4_3*m_4_3>. This is done on condor
3. hadd all the output rootfiles 
4. Run corrections (CBWC and PileUp Corr.)
5. Calculate the systematic uncertainty 
6. Plot the results


## How to run the code ##

You will need to first submit the jobs by going into `scripts` and running:

` submit_full_analysis.sh `

Wait for the jobs to finish (~1-3 hours).

Then run the script in ` scripts ` that adds the histograms, applies the pileup and CWBC correction and plots the result:

` ./run_corr_and_plot.sh `

an image of the plot cumulants will show up in `imgs`. 




 
 
