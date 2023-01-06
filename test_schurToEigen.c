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
    double norm = hqr(n,n,0,n - 1,T,eigValsReal,eigValsImag,1,eigenMat);
    if (norm < 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %d\n",-norm);
        return 1;
    }

    schurToEigen(0, n - 1, norm, n, eigValsReal, eigValsImag, T, eigenMat);

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
        // Check if the current eigenvalue is real or imaginary.
        // If it is real, we need only check Ax = ax
        if (eigValsImag[k] == 0.0) {
            // Real eigenvalue (eigenValsReal[k]) with eigenvector eigenMat(:,k)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Ax[i] += a0(i,j) * eigenMat0(j,k);
                }
                ax[i] = eigValsReal[k] * eigenMat0(i,k);
            }
        } else if (eigValsImag[k] > 0.0) {
            // Imaginary with positive imaginary part so has associated eigenvector
            // eigenMat(:,k) + i * eigenMat(:,k+1)
            // compute Ax,Ay
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Ax[i] += a0(i,j) * eigenMat0(j,k);
                    Ay[i] += a0(i,j) * eigenMat0(j, k + 1);
                }
                ax[i] = eigValsReal[k] * eigenMat0(i,k) - eigValsImag[k] * eigenMat0(i,k+1);
                ay[i] = eigValsReal[k] * eigenMat0(i,k + 1) + eigValsImag[k] * eigenMat0(i,k);
            }
        } else {
            // Imaginary with negative imaginary part so has associated eigenvector
            // eigenMat(:,k-1) - i * eigenMat(:,k)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Ax[i] += a0(i,j) * eigenMat0(j,k - 1);
                    Ay[i] -= a0(i,j) * eigenMat0(j, k);
                }
                ax[i] = eigValsReal[k] * eigenMat0(i, k - 1) + eigValsImag[k] * eigenMat0(i, k);
                ay[i] = eigValsImag[k] * eigenMat0(i, k - 1) - eigValsReal[k] * eigenMat0(i, k);
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
            normR += fabs(diffMatReal[ i + j * n]) + fabs(diffMatImag[i + j * n]);
            normA += fabs(a0(i,j));
            printf("%1.10e ",diffMatImag[i + j * n]);
        }
        printf("\n");
    }
    printf("%1.5e\n%1.5e\n",normR,normA);
    printf("n = %6d, Eigenvector Check: %1.10e\n", n, normR / normA);

    free(A);
    free(T);
    free(eigenMat);
    return 0;
}
