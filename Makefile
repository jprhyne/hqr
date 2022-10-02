CC = gcc
FC = gfortran

all: main_hqr_test.exe main_hqr_goto.exe

hqr.o: hqr.f
	$(FC) -c hqr.f

main_hqr_test.o: main_hqr_test.c
	$(CC) -c main_hqr_test.c -lm

main_hqr_goto.o: main_hqr_goto.c
	$(CC) -c main_hqr_goto.c -lm

main_hqr_goto.exe: main_hqr_goto.o hqr.o
	$(FC) hqr.o main_hqr_goto.o -o $@

main_hqr_test.exe: main_hqr_test.o hqr.o
	$(FC) hqr.o main_hqr_test.o -o $@

clean:
	rm -f *.exe *.o *.out

