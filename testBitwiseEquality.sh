#!/bin/env bash
# Call make as it only does something if the file updates
make "testBitwiseEquality.exe"

if [ -f "bitEq.txt" ]; then
	rm "bitEq.txt"
fi

for n in {10..1010..100}; do
	./testBitwiseEquality.exe -n $n
done
