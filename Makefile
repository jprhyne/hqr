include make.inc

all: test_schurVectors.exe test_hqr2schur_fortran.exe test_hqr2eigen_fortran.exe testBitwiseEquality.exe test_schurToEigen.exe

cdivF.o: cdiv.f
	$(FC) -c cdiv.f -g -o cdivF.o

cdivC.o: cdiv.c
	$(CC) -c $^ -o cdivC.o -g

doubleSubDiagonalSearch.o: doubleSubDiagonalSearch.c 
	$(CC) -c $^ -g

formShift.o: formShift.c
	$(CC) -c $^ -g

hqr2Schur.o: hqr2Schur.f
	$(FC) -c $^ -g

hqr2Eigen.o: hqr2Eigen.f
	$(FC) -c $^ -g

hqrC.o: hqr.c
	$(CC) -c hqr.c -g -o hqrC.o

hqr.o: hqr.f
	$(FC) -c hqr.f -g

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

schurToEigen.o: schurToEigen.c
	$(CC) -c $^ -g

testBitwiseEquality.o: testBitwiseEquality.c
	$(CC) -c $^ -lm -g

test_hqr2schur_fortran.o: test_hqr2schur_fortran.c
	$(CC) -c $^ -lm -g

test_hqr2eigen_fortran.o: test_hqr2eigen_fortran.c
	$(CC) -c $^ -lm -g

test_schurVectors.o: test_schurVectors.c
	$(CC) -c $^ -lm -g

test_schurToEigen.o: test_schurToEigen.c
	$(CC) -c $^ -lm -g

test_schurVectors.exe: test_schurVectors.o formShift.o hqrC.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o cdivC.o qrIterationVec.o matmul.o matsub.o
	$(LOADER) -o $@ $^

test_schurToEigen.exe: test_schurToEigen.o formShift.o hqrC.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o cdivC.o qrIterationVec.o matmul.o matsub.o schurToEigen.o
	$(LOADER) -o $@ $^

test_hqr2schur_fortran.exe: test_hqr2schur_fortran.o hqr2Schur.o hqr2Eigen.o cdivC.o matmul.o matsub.o cdivF.o
	$(LOADER) -o $@ $^

test_hqr2eigen_fortran.exe: test_hqr2eigen_fortran.o hqr2Schur.o hqr2Eigen.o cdivC.o matmul.o matsub.o cdivF.o formShift.o hqrC.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o cdivC.o qrIterationVec.o
	$(LOADER) -o $@ $^

testBitwiseEquality.exe: testBitwiseEquality.o formShift.o hqrC.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o cdivC.o qrIterationVec.o hqr2Schur.o hqr2Eigen.o cdivF.o schurToEigen.o
	$(LOADER) -o $@ $^

clean:
	rm -f *.exe *.o *.out .*swp

