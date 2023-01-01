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
    printf("test_hqr2eigen_fortran.exe [-n sizeOfMatrix | -v | -h]\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-v: verbose flag that prints out the results of eispack hqr\n");
    printf("\t\tBy default, nothing is printed to the console\n");
    printf("\t-t: testing flag that only prints the expected vs the");
    printf("\t\tactual computed eigenvalues.");
    printf("\t-h: Print this help dialogue\n");
    printf("\t--jobv: Flag that tells us to compute the eigenvectors");
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
		int testFlag = 0;
        int eigenVectorFlag = 0;

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
	    } else if ( strcmp ( *(argv + i), "--jobv") == 0) {
                eigenVectorFlag=1;
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
    // Getting here means that we have successfully ran all of 
    // hqr and got an answer. 
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

    //  check || A * Z - Z * T  ||_F / || A ||_F
    double normR, normA;
    normR = 0.0e+00;
    // Yes, this is lazy and inefficient, however it 
    // accomplishes our goals in a reasonable time frame
    double *lhs = matmul(A,n,n,Z,n,n);
    double *rhs = matmul(Z,n,n,T,n,n);
    double *ans = matsub(lhs,n,n,rhs,n,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            normR += ans[i + j * n] * ans[i + j * n];
    normR = sqrt( normR );
    normA = 0.0e+00;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            normA += a0(i,j) * a0(i,j) ;
        }
    }
    normA = sqrt( normA );
    
    // Now we check if our eigen matrix is invertible
    // Since finding the inverse is expensive, and inaccurate, we will get around
    // The invertability check by finding the eigenvalues and checking to make 
    // sure none of them are 0, then we check representation.
    //
    // This assumes that hqr works correctly
    double *eigEigValReal = malloc(n * sizeof(double));
    double *eigEigValImag = malloc(n * sizeof(double));
    hqr(n,n,1,n,Z,eigEigValReal,eigEigValImag,0,NULL);
    // Check if any eigenvalue is 0
    for (int k = 0; k < n; k++) 
        if (eigEigValReal[k] == 0 && eigEigValImag[k] == 0)
            printf("eigenValue of the Zrix is 0 at index %4d",k);
    printf("%% [ REPRES ] n = %4d; check = [ %1.10e ];\n", n, normR / normA );
    /*
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%1.10e",t0(i,j));
        }
        printf("\n");
    }
    */
    free(A);
    free(T);
    free(Z);
    free(rhs);
    free(lhs);
    free(ans);
    return 0;
}
