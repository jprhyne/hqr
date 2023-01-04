#!/bin/env bash
# This script requires either a linux based operating system
# or a Unix one with bash installed. May need to change the shebang to get this
# to work

# These following values are the dimensions we want to test. If
# the user passes N as the dimension we test 1,...,N dimensions
# unless the user also specifies the number of steps
numIter=50
numStep=1

usage () {
	echo "test_hqr_conversion.sh [-h|-n #|-s #]"
	echo "-n is followed by an integer. This determines the maximum dimension we are testing."
	echo "-s is followed by an integer. This determines the space in between dimensions we test."
	echo "   If this is greater than the number of steps, then it is reassigned to the number of steps"
	exit 0
}
while test $# -gt 0; do
	case "$1" in 
		-h)
			usage
			;;
		-n)
			shift
			tmp=$1
			if [  ${tmp##*[!0-9]} ]; then
				# tmp is an integer so replace the numIter variable 
				numIter=$tmp
			fi
			;;
		-s) 
			shift
			tmp=$1
			if [ ${tmp##*[!0-9]} ]; then
				# tmp is an integer, so replace the numStep variable
				numStep=$tmp
			fi
			;;
	esac
	shift
done
if [ "$numStep" -gt "$numIter" ]; then
	numStep=$numIter
fi

# If the executable main_hqr_loopByLoopConversion.exe is not present, then we run make to compile it
if [ ! -f "main_hqr_loopByLoopConversion.exe" ]; then
	make
fi
# Check if the file "results.txt" exists. If so, ask if the user
# wants to overwrite it
if [ -f "results.txt" ]; then
	read -p "results.txt exists. Do you wish to overwrite? [y/N] " response
	if [ "$response" == "y" ] || [ "$response" == "Y" ]; then
		rm results.txt
	else 
		echo "Move the file or rename it."
		exit 1
	fi
fi
echo "" > results.txt
for k in $(seq 1 "$numStep" "$numIter"); do
	echo -n "k=$k, " >> results.txt
	./main_hqr_loopByLoopConversion.exe -t -n $k >> results.txt
done
