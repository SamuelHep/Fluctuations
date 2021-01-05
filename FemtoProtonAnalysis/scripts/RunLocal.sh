#!/bin/bash

FILELIST=../input_parameters/smallSample.list
parameterList=../input_parameters/parameter.list
outfilename=../rootfiles/cumulants/smallSample_out.root

root -b -l -q ../macros/RunEventsLocal.C\(\"$FILELIST\",\"$parameterList\",\"$outfilename\"\)
