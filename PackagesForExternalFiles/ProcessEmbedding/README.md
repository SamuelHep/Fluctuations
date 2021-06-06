# Proton Embedding #
If you'd like to generate your own efficiency files, use this package to make a reduced 
proton embedding file.

Simply run the command `./ProcessEmbedding.sh` and outfiles will be generated in 
`$output_parent_directory/embed/out`

After finishing, go to `FemtoProtonAnalysis/scripts/` and run `./CalcTPCEfficiency.sh`

For TOF matching efficiency, go to `FemtoProtonAnalysis/scripts/` and run `./CalcTOFEfficiency.sh`

Both scripts will generate efficiency files in `FemtoProtonAnalysis/rootfiles/eff/`