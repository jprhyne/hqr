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
			if 
			fi
			numIter=$1
			;;
		-s) 
			shift
			numIter=$1
			;;
	esac
	shift
done
echo "numIter=$numIter"
echo "numStep=$numStep"
