include make.inc

all: main_hqr_loopByLoopConversion.exe test_schurVectors.exe test_hqr2_fortran.exe testBitwiseEquality.exe

cdivF.o: cdiv.f
	$(FC) -c cdiv.f -g -o cdivF.o

cdivC.o: cdiv.c
	$(CC) -c $^ -o cdivC.o -g

doubleSubDiagonalSearch.o: doubleSubDiagonalSearch.c 
	$(CC) -c $^ -g

formShift.o: formShift.c
	$(CC) -c $^ -g

hqr2.o: hqr2.f
	$(FC) -c hqr2.f -g

hqrC.o: hqr.c
	$(CC) -c hqr.c -g -o hqrC.o

hqr.o: hqr.f
	$(FC) -c hqr.f -g

main_hqr_loopByLoopConversion.o: main_hqr_loopByLoopConversion.c
	$(CC) -c main_hqr_loopByLoopConversion.c -lm -g

matmul.o: matmul.c
	$(CC) -c $^ -g

matsub.o: matsub.c
	$(CC) -c $^ -g

qrIteration.o: qrIteration.c
	$(CC) -c $^ -g -lm

qrIterationVec.o: qrIterationVec.c
	$(CC) -c $^ -g -lm

subDiagonalSearch.o: subDiagonalSearch.c
	$(CC) -c $^ -g

testBitwiseEquality.o: testBitwiseEquality.c
	$(CC) -c $^ -lm -g

test_hqr2_fortran.o: test_hqr2_fortran.c
	$(CC) -c $^ -lm -g

test_schurVectors.o: test_schurVectors.c
	$(CC) -c $^ -lm -g

main_hqr_loopByLoopConversion.exe: main_hqr_loopByLoopConversion.o formShift.o hqr.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o hqr2.o cdivF.o cdivC.o qrIterationVec.o
	$(LOADER) -o $@ $^

test_schurVectors.exe: test_schurVectors.o formShift.o hqrC.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o cdivC.o qrIterationVec.o matmul.o matsub.o
	$(LOADER) -o $@ $^

test_hqr2_fortran.exe: test_hqr2_fortran.o hqr2.o cdiv.o matmul.o matsub.o
	$(LOADER) -o $@ $^

testBitwiseEquality.exe: testBitwiseEquality.o formShift.o hqrC.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o cdivC.o qrIterationVec.o hqr2.o
	$(LOADER) -o $@ $^

clean:
	rm -f *.exe *.o *.out .*swp

