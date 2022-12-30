#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define a0(i,j) A[(i) + (j) * n]
#define b0(i,j) B[(i) + (j) * n]
#define z0(i,j) z[(i) + (j) * n]

// This may be bad practice, but we put all malloc'd entities
// as global variables in order to make freeing the
// memory easier
double* A;
double* B;
double* wr;
double* wi;
double* z;
double* eigenValsReal;
double* eigenValsImag;

extern double*matmul(double*A, int nA, int mA, double *B, int nB, int mB);

extern double*matsub(double*A, int nA, int mA, double *B, int nB, int mB);

extern void hqr2_(int *nm,int *n,int *low,int *igh, double *h, double *wr,
        double *wi, double *z, int *ierr);

void usage()
{
    printf("test_hqr2_fortran.exe [-n sizeOfMatrix | -v | -h]\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-v: verbose flag that prints out the results of eispack hqr\n");
    printf("\t\tBy default, nothing is printed to the console\n");
    printf("\t-t: testing flag that only prints the expected vs the");
    printf("\t\tactual computed eigenvalues.");
    printf("\t-h: Print this help dialogue\n");
    printf("\t--jobv: Flag that tells us to compute the eigenvectors");
}

void freeMemory()
{
    free(A);
    free(B);
    free(wr);
    free(wi);
    free(z);
    free(eigenValsReal);
    free(eigenValsImag);
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
        // arrays that will hold the differences in the eigenvalues of hqr
        // and this implementation
	double eigRealDiff[n];
	double eigImagDiff[n];
        double eigVec[n*n];

	// Allocate the memory for A to be generated. It will contain n^2 
	// elements where each element is a double precision floating point number
	A = (double *) malloc( n * n *  sizeof(double));
	B = (double *) malloc( n * n *  sizeof(double));
	// Create a vector to store the real parts of the eigenvalues
	wr = (double *) malloc( n *  sizeof(double));
	// Create a vector to store the real parts of the eigenvalues
	wi = (double *) malloc( n *  sizeof(double));

        z = (double *) malloc( n * n * sizeof(double));

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
            }
        }

        // Store a copy of A into B so that we can run our own
        // version of hqr written in C 
        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                b0(i,j) = a0(i,j);
	hqr2_( &n, &n, &ione, &n, A, wr, wi, z, &ierr);
	/*
	 * Below prints out wr and wi for inspection via 
	 * MATLAB/visual inspection
	 */
	// Check that V^T V = I 
    double orthZ, tmp;
    printf("%% [ ORTH ] n = %4d; checks = [ ", n );
    orthZ = 0e+00;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
            for (int k = 0; k < n; k++) {
                tmp -= z[k + i *n]*z[k + j * n];
            }
            orthZ += tmp * tmp;
        }
    }
    orthZ = sqrt( orthZ );
    printf(" %1.10e", orthZ );

	// zero-out below quasi diagonal of H by looking wi
    // First, zero out everything below the 1st subdiagonal
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < i - 1; j++) 
            A[i + j * n] = 0;
    // if eigValsImag[k]  = 0 then the sub diagonal elements need to be 0
    // If eigValsImag[k] != 0 then we have a schur block  
    int k;
    for (k = 0; k < n-1; k++) {
        if (wi[k] == 0) {
            A[(k + 1) + (k) * n] = 0;
        } else if (k < n-2){
            // This means we are in a schur block, so the next sub diagonal
            // element must be 0
            A[(k + 2) + (k + 1) * n] = 0;
            k++;
        }
    }
	// Check that A = V * H * V^T
    double normR, normA;
    normR = 0.0e+00;
    double *Zt = malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Zt[i + j * n] = z[j + i * n];
    // Yes, this is lazy and inefficient, however it 
    // accomplishes our goals in a reasonable time frame
    double *ZT = matmul(z,n,n,A,n,n);
    double *rhs = matmul(ZT,n,n,Zt,n,n);
    double *ans = matsub(rhs,n,n,B,n,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            normR += ans[i + j * n] * ans[i + j * n];
    normR = sqrt( normR );
    normA = 0.0e+00;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            normA += B[i + j * n] * B[i + j * n];
        }
    }
    normA = sqrt( normA );
    printf(" %1.10e", normR / normA );
    printf(" ];\n");
    free(A);
    free(B);
    free(Zt);
    free(ZT);
    free(rhs);
    free(ans);
    return 0;
}
