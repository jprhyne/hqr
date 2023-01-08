#!/bin/env bash
numIter=50
numStep=1

for k in $(seq 1 "$numStep" "$numIter"); do
	for seed in $(seq 1 "$numStep" "$numIter"); do
		echo -n "seed=$seed,   " 
		./test_hqr2eigen_fortran.exe -n $k -s $seed
	done
done
