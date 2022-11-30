#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define a0(i,j) A[(i) + (j) * n]
#define b0(i,j) A[(i) + (j) * n]

extern int hqr(int nm, int n, int low, int igh, double *A, double *eigenValsReal, double *eigenValsImag, int schurVectorFlag);

extern double*matmul(double*A, int nA, int mA, double *B, int nB, int mB);

extern double*matsub(double*A, int nA, int mA, double *B, int nB, int mB);

void usage()
{
    printf("main_hqr_test.exe [-h | -n sizeOfMatrix | -s seed]\n");
    printf("\t-h: Print this help dialogue\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-s: Sets the seed. The following argument must be a positive integer.\n");
    printf("\t\tThe default value is 28.");
}

int main (int argc, char **argv) 
{
    // Seeds the random number generator for repeatability
    int seed = 28;
    // Default size of the matrix A
    int n = 20;
    for(int i = 1; i < argc; i++){
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
        } else if ( strcmp ( *(argv + i), "-s") == 0) {
            seed = atoi( *(argv + i + 1) ); 
            if (seed <= 0)
                usage();
            i++;
        } 
    }
    srand(seed);
    // First create a matrix A as upper hessenberg matrix.
    double *A = calloc(n*n,sizeof(double));
    // These values are not needed for this particular file's
    // testing, however this is needed to run hqr.c
    double *eigValsReal = malloc(n*n*sizeof(double));
    double *eigValsImag = malloc(n*n*sizeof(double));
	// Generate A as a random matrix.
    for(int i = 0; i < n; i++) {
        int start = 0;
        if (i - 1 > 0)
            start = i - 1;
 	    for(int j = start; j < n; j++) {
            double val = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	        a0(i,j) = val; 
        }
    }
    // Now we call hqr. At the end A will contain the schur vectors
    int ret = hqr(n,n,1,n,A,eigValsReal,eigValsImag,1);
    if (ret != n) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %d\n",ret);
        return 1;
    }
    // Getting here means that we have successfully ran all of 
    // hqr and got an answer, so now we check if our Schur vectors are correct

    //First, we have to compute A^\top
    double *B = calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++)
            b0(j,i) = a0(i,j);
    double* C = matmul(A,n,n,B,n,n);
    // Make an identity matrix
    double* eye = calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i ++) 
        eye[i + i * n] = 1; 
    double* diffMat = matsub(C,n,n,eye,n,n);
    // compute the sum of the absolute value of the elements of diffMat
    double norm = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norm += diffMat[i + j*n];
    printf("Sum of the absolute elements of diffmat: %1.20f",norm);
}
