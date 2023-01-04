#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include"externalFunctions.h"
extern void schurToEigen(int low, int igh, double norm, int n, double *eigenValsReal,
            double *eigenValsImag, double *T, double *eigenMatrix);
#define a0(i,j) A[(i) + (j) * n]
#define eigenMat0(i,j) eigenMat[(i) + (j) * n]
#define t0(i,j) T[(i) + (j) * n]

void usage()
{
    printf("test_schurToEigen.exe [-h | -n sizeOfMatrix | -s seed]\n");
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
    double *A = (double *) calloc(n*n,sizeof(double));
    double *T = (double *) calloc(n*n,sizeof(double));
    // These values are not needed for this particular file's
    // testing, however this is needed to run hqr.c
    double *eigValsReal = (double *) malloc(n*n*sizeof(double));
    double *eigValsImag = (double *) malloc(n*n*sizeof(double));
	// Generate A as a random matrix.
    for(int i = 0; i < n; i++) {
        int start = 0;
        if (i - 1 > 0)
            start = i - 1;
 	    for(int j = start; j < n; j++) {
            double val = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	        a0(i,j) = val; 
	        t0(i,j) = val; 
        }
    }
    // Create a matrix to store the Schur Vectors
    double *eigenMat = (double *) calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i++)
        eigenMat[i + i * n] = 1;
    // Now we call hqr. At the end eigenMat will contain the schur vectors
    double norm = hqr(n,n,1,n,T,eigValsReal,eigValsImag,1,eigenMat);
    if (norm < 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %d\n",-norm);
        return 1;
    }

    schurToEigen(1, n, norm, n, eigValsReal, eigValsImag, T, eigenMat);

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
    double *lhs = matmul(A,n,n,eigenMat,n,n);
    double *rhs = matmul(eigenMat,n,n,T,n,n);
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
    
    // Now we check if our eigenMatrix is invertible
    // Since finding the inverse is expensive, and inaccurate, we will get around
    // The invertibility check by finding the eigenvalues and checking to make 
    // sure none of them are 0, then we check representation.
    //
    // This assumes that hqr works correctly
    double *eigEigValReal = malloc(n * sizeof(double));
    double *eigEigValImag = malloc(n * sizeof(double));
    hqr(n,n,1,n,eigenMat,eigEigValReal,eigEigValImag,0,NULL);
    // Check if any eigenvalue is 0
    for (int k = 0; k < n; k++) 
        if (eigEigValReal[k] == 0 && eigEigValImag[k] == 0)
            printf("eigenValue of the eigenMatrix is 0 at index %4d",k);
    printf("%% [ REPRES ] n = %4d; check = [ %1.10e ];\n", n, normR / normA );
    free(A);
    //free(T);
    free(eigenMat);
    free(rhs);
    free(lhs);
    free(ans);
    return 0;
}
