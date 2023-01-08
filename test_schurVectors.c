#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include "externalFunctions.h"
#define a0(i,j) A[(i) + (j) * n]
#define schurMat0(i,j) schurMat[(i) + (j) * n]
#define t0(i,j) T[(i) + (j) * n]

void usage()
{
    printf("test_schurVectors.exe [-h | -n sizeOfMatrix | -s seed | -t]\n");
    printf("\t-h: Print this help dialogue\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-s: Sets the seed. The following argument must be a positive integer.\n");
    printf("\t\tThe default value is 28.");
    printf("\t-t: This flag tells us if we want the output in a human readable format\n");
    printf("\t\twe default to machine readable to allow for easier plot creation");
}

int main (int argc, char **argv) 
{
    // Seeds the random number generator for repeatability
    int seed = 28;
    // Default size of the matrix A
    int n = 20;
    //Flag that determines if we want the output in a 
    //human readable format or in machine readable format
    //We default to machine
    int testFlag = 0;
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
        } else if ( strcmp( argv[i], "-t") == 0) {
            testFlag = 1;
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
    double *schurMat = (double *) calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i++)
        schurMat[i + i * n] = 1;
    // Now we call hqr. At the end schurMat will contain the schur vectors
    double norm = hqr(n,n,0,n-1,T,eigValsReal,eigValsImag,1,schurMat);
    if (norm < 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %e\n",-norm);
        return 1;
    }
    /*
    for (int i = 0; i < n; i++){
        for(int j=0;j<n;j++) {
            printf("%1.11f", schurMat[i + j * n]);
        }
        printf("\n");
    }
    */
    // Getting here means that we have successfully ran all of 
    // hqr and got an answer, so now we check if our Schur vectors are correct
    //  check || Z' * Z - I ||_F
    double orthZ, tmp;
    orthZ = 0e+00;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
            for (int k = 0; k < n; k++) {
                tmp -= schurMat0(k,i)*schurMat0(k,j);
            }
            orthZ += tmp * tmp;
        }
    }
    orthZ = sqrt( orthZ );

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

    //  check || A - Z * T * Z^T ||_F / || A ||_F
    double normR, normA;
    normR = 0.0e+00;
    double *Zt = malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Zt[i + j * n] = schurMat0(j,i);
    // Yes, this is lazy and inefficient, however it 
    // accomplishes our goals in a reasonable time frame
    double *ZT = matmul(schurMat,n,n,T,n,n);
    double *rhs = matmul(ZT,n,n,Zt,n,n);
    double *ans = matsub(rhs,n,n,A,n,n);
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
    if (testFlag)
        printf("%% [ hqr2schur C ] n = %8d; seed = %8d; checks = [ %8.2e %8.2e ];\n", n, seed, orthZ, normR/normA);
    else
        printf( "%8d %8d %6.1e %6.1e\n", n, seed, orthZ, normR/normA);
    free(A);
    free(T);
    free(schurMat);
    free(Zt);
    free(ZT);
    free(rhs);
    free(ans);
    return 0;
}
