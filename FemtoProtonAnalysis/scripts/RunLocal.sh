#!/bin/bash

FILELIST=../input_parameters/femto_3GeV_test.list
parameterList=../input_parameters/parameter_n0p5_0_norm.list
outfilename=test_bootstrap.root

root -b -l -q ../macros/RunEventsLocal.C\(\"$FILELIST\",\"$parameterList\",\"$outfilename\"\)
