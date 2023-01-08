#!/bin/env bash
numIterN=100
numStepN=1

numIterSeed=50
numStepSeed=1

for k in $(seq 1 "$numStepN" "$numIterN"); do
	for seed in $(seq 1 "$numStepSeed" "$numIterSeed"); do
		echo -n "seed=$seed,   " 
		#./test_hqr2eigen_fortran.exe -n $k -s $seed
		./test_schurVectors.exe -n $k -s $seed
	done
done
