CC = gcc
FC = gfortran

all: main_hqr_loopByLoopConversion.exe

hqr.o: hqr.f
	$(FC) -c hqr.f -g

hqr2.o: hqr2.f
	$(FC) -c hqr2.f -g

cdiv.o: cdiv.f
	$(FC) -c cdiv.f -g

hqr_destructive.o: hqr_destructive.f
	$(FC) -c hqr_destructive.f -g

formShift.o: formShift.c
	$(CC) -c $^ -g

subDiagonalSearch.o: subDiagonalSearch.c
	$(CC) -c $^ -g

doubleSubDiagonalSearch.o: doubleSubDiagonalSearch.c 
	$(CC) -c $^ -g

qrIteration.o: qrIteration.c
	$(CC) -c $^ -g -lm

main_hqr_test.o: main_hqr_test.c
	$(CC) -c main_hqr_test.c -lm -g

main_hqr_goto.o: main_hqr_goto.c
	$(CC) -c main_hqr_goto.c -lm -g

main_hqr_goto.exe: main_hqr_goto.o hqr.o
	$(FC) hqr.o main_hqr_goto.o -o $@

main_hqr_loopByLoopConversion.o: main_hqr_loopByLoopConversion.c
	$(CC) -c main_hqr_loopByLoopConversion.c -lm -g

main_hqr_loopByLoopConversion.exe: main_hqr_loopByLoopConversion.o formShift.o hqr.o subDiagonalSearch.o doubleSubDiagonalSearch.o qrIteration.o hqr2.o cdiv.o
	$(FC) -o $@ $^

main_hqr_test.exe: main_hqr_test.o hqr.o
	$(FC) hqr.o main_hqr_test.o -o $@

clean:
	rm -f *.exe *.o *.out .*swp

