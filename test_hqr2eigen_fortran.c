#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "externalFunctions.h"
#define a0(i,j) A[(i) + (j) * n]
#define t0(i,j) T[(i) + (j) * n]
#define z0(i,j) Z[(i) + (j) * n]

// This may be bad practice, but we put all malloc'd entities
// as global variables in order to make freeing the
// memory easier
double* A;
double* T;
double* Z;
double* eigValsReal;
double* eigValsImag;

void usage()
{
    printf("test_hqr2eigen_fortran.exe [-h | -n sizeOfMatrix | -s seed | -t]\n");
    printf("\t-h: Print this help dialogue\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-s: Sets the seed. The following argument must be a positive integer.\n");
    printf("\t\tThe default value is 28.");
    printf("\t-t: This flag tells us if we want the output in a human readable format\n");
    printf("\t\twe default to machine readable to allow for easier plot creation");
}

int main(int argc, char ** argv) {
/*
 * Declaring variables to pass into the hqr_ subroutine
 */
        int i, ierr, j, n;

        // Store the constant 1 as a variable to be sent into hqr_
        int ione=1;
        
        /*
         * flag that determines if we want to print the results
         * of hqr_ to the console
         * By default, we do not do this
         */
        int printFlag = 0;
        // This prints results in a human readable way
        int testFlag = 0;


        // Seeds the random number generator for repeatability
        // Old seed was 734
        int seed = 28;

        // Default size of the matrix A
        n = 20;

        // Checking for if the user wants to set the size of the matrix A
        for(i = 1; i < argc; i++){
            char *argument = argv[i];
            if (strcmp(argument, "-h") == 0) {
                usage();
                return 0;
            } else if( strcmp( *(argv + i), "-n") == 0) {
                // Grab the next term as this is the size
                n  = atoi( *(argv + i + 1) );
                //If the size is non-positive, exit immediately
                if ( n <= 0 )
                    usage();
                // Increment i to skip over this number
                i++;
            } else if ( strcmp ( *(argv + i), "-v") == 0) {
                printFlag=1;
            } else if ( strcmp ( *(argv + i), "-t") == 0) {
                testFlag=1;
            } else if ( strcmp ( *(argv + i), "-s") == 0) {
                seed = atoi( *(argv + i + 1) ); 
                if (seed <= 0)
                    usage();
                i++;
            }
        }
        //Uncomment if we want to test the output of qrIteration.c
        //testingFile = fopen("outputFileC.txt","w");
        srand(seed);

        // Allocate the memory for A to be generated. It will contain n^2 
        // elements where each element is a double precision floating point number
        A = (double *) malloc( n * n *  sizeof(double));
        T = (double *) malloc( n * n *  sizeof(double));
        // Create a vector to store the real parts of the eigenvalues
        eigValsReal = (double *) malloc( n *  sizeof(double));
        // Create a vector to store the real parts of the eigenvalues
        eigValsImag = (double *) malloc( n *  sizeof(double));

        Z = (double *) malloc( n * n * sizeof(double));

        for (i = 0; i < n; i++) 
            z0(i,i) = 1;

        // Generate A as a random matrix.
        for(i = 0; i < n; i++) {
            int start = 0;
            if (i - 1 > 0)
                start = i - 1;
        for(j = start; j < n; j++) {
                double val = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
                a0(i,j) = val; 
                t0(i,j) = val;
            }
        }
        hqr2eigen_( &n, &n, &ione, &n, T, eigValsReal, eigValsImag, Z, &ierr);
        // Zero out below quasi diagonal elements of T
        // First, zero out everything below the 1st subdiagonal
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < i - 1; j++) 
                t0(i,j) = 0;
        // if eigValsImag[k]  = 0 then the sub diagonal elements need to be 0
        // If eigValsImag[k] != 0 then we have a schur block  
        int k;
        for (k = 0; k < n-1; k++) {
            if (eigValsImag[k] == 0) {
                t0(k+1,k) = 0;
            } else if (k < n-2){
                // This means we are in a schur block, so the next sub diagonal
                // element must be 0
                t0(k+2,k+1) = 0;
                k++;
            }
        }
        double diffMatReal[n * n];
        double diffMatImag[n * n];
        for (k = 0; k < n; k++) {
            // Redeclare on each loop to not have to reset each value to 0
            // Vectors that will store Av = Ax + i * Ay
            double Ax[n];
            double Ay[n];
            // Vectors that will store \lambda v = (a + ib)(x+iy) = (ax - by) + i(ay + bx)
            double ax[n];
            double ay[n];
            for (int l = 0; l < n; l++) {
                Ax[l] = 0.0;
                Ay[l] = 0.0;
                ax[l] = 0.0;
                ay[l] = 0.0;
            }
            // Check if the current eigenvalue is real or imaginary.
            // If it is real, we need only check Ax = ax
            if (eigValsImag[k] == 0.0) {
                // Real eigenvalue (eigenValsReal[k]) with eigenvector eigenMat(:,k)
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        Ax[i] += a0(i,j) * z0(j,k);
                    }
                    ax[i] = eigValsReal[k] * z0(i,k);
                }
            } else if (eigValsImag[k] > 0.0) {
                //continue;
                // Imaginary with positive imaginary part so has associated eigenvector
                // eigenMat(:,k) + i * eigenMat(:,k+1)
                // compute Ax,Ay
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        Ax[i] += a0(i,j) * z0(j,k);
                        Ay[i] += a0(i,j) * z0(j, k + 1);
                    }
                    ax[i] = eigValsReal[k] * z0(i,k) - eigValsImag[k] * z0(i,k+1);
                    ay[i] = eigValsReal[k] * z0(i,k + 1) + eigValsImag[k] * z0(i,k);
                }
            } else {
                //continue;
                // Imaginary with negative imaginary part so has associated eigenvector
                // eigenMat(:,k-1) - i * eigenMat(:,k)
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        Ax[i] += a0(i,j) * z0(j,k - 1);
                        Ay[i] -= a0(i,j) * z0(j, k);
                    }
                    ax[i] = eigValsReal[k] * z0(i, k - 1) + eigValsImag[k] * z0(i, k);
                    ay[i] = eigValsImag[k] * z0(i, k - 1) - eigValsReal[k] * z0(i, k);
                }
            } 
            // Now we construct the kth column of diffMatReal and diffMatImag. Note that if we had a real eigenvalue the
            // kth column of diffMatImag will be exactly 0
            for (int i = 0; i < n; i++) {
                // We want the kth column of diffMatReal to be Ax - ax
                // and the kth column of diffMatImag to be Ay - ay
                diffMatReal[i + k * n] = Ax[i] - ax[i];
                diffMatImag[i + k * n] = Ay[i] - ay[i];
            }
        }
        // Now we have computed AV - VD and stored the real and imaginary parts separately
        // Since we want to get a relative error we must compute \|AV - VD\| in order
        // for this to be as close to accurate as possible we will compute the sum 
        // of the absolute elements of AV - VD using the fact that |a + bi| = sqrt(a^2 + b^2)
        double normR = 0.0;
        double normA = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++){
                double tmpA = diffMatReal[i + j * n];
                double tmpB = diffMatImag[i + j * n];
                normR += sqrt(tmpA * tmpA + tmpB * tmpB); 
                normA += fabs(a0(i,j));
            }
        }
        if (testFlag)
            printf("n = %6d, Eigenvector Check: %1.10e\n", n, normR / normA);
        else
            printf( "%d %d %1.10e\n", n, seed, normR/normA);
        free(A);
        free(T);
        free(Z);
        return 0;
}
