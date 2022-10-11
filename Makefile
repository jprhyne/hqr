CC = gcc
FC = gfortran

all: main_hqr_loopByLoopConversion.exe

hqr.o: hqr.f
	$(FC) -c hqr.f -g

hqr_destructive.o: hqr_destructive.f
	$(FC) -c hqr_destructive.f -g

hqr_formshift.o: hqr_formshift.f
	$(FC) -c hqr_formshift.f -g

main_hqr_test.o: main_hqr_test.c
	$(CC) -c main_hqr_test.c -lm -g

main_hqr_goto.o: main_hqr_goto.c
	$(CC) -c main_hqr_goto.c -lm -g

main_hqr_goto.exe: main_hqr_goto.o hqr.o
	$(FC) hqr.o main_hqr_goto.o -o $@

main_hqr_loopByLoopConversion.o: main_hqr_loopByLoopConversion.c
	$(CC) -c main_hqr_loopByLoopConversion.c -lm -g

main_hqr_loopByLoopConversion.exe: main_hqr_loopByLoopConversion.o hqr_formshift.o hqr.o
	$(FC) hqr_formshift.o main_hqr_loopByLoopConversion.o hqr.o -o $@

main_hqr_test.exe: main_hqr_test.o hqr.o
	$(FC) hqr.o main_hqr_test.o -o $@

clean:
	rm -f *.exe *.o *.out

