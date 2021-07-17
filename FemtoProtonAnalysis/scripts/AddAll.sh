#!/bin/bash

sysArray=("norm" "SYS1" "SYS2" "SYS3" "SYS4" "SYS5" "SYS6" "SYS7" "SYS8")
rapArray=("n0p2_0" "n0p3_0" "n0p4_0" "n0p5_0")
ptArray=("n0p5_0_pt1" "n0p5_0_pt2" "n0p5_0_pt3" "n0p5_0_pt4") 

for sys in ${sysArray[@]}; do
    for rap in ${rapArray[@]}; do
	python MadAdder.py ${rap}_${sys} profiles_${rap}_${sys}.root
	wait
    done 
    for pt in ${ptArray[@]}; do
	python MadAdder.py ${pt}_${sys} profiles_${pt}_${sys}.root
	wait
    done 

done 
