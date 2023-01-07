#!/bin/env bash
# Call make as it only does something if the file updates
make "testBitwiseEquality.exe"

if [ -f "bitEq.txt" ]; then
	rm "bitEq.txt"
fi

for n in {10..1010..10}; do
	# The & allows for the OS to parallelize these jobs. Could cause issues with writing to the file
	# however this is not an issue as long as your implementation of C's fprintf is POSIX compliant
	./testBitwiseEquality.exe -n $n &
done
