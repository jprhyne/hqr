#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#define a0(i,j) A[(i) + (j) * n]
#define schurMat0(i,j) schurMat[(i) + (j) * n]
#define t0(i,j) T[(i) + (j) * n]

extern int hqr(int nm, int n, int low, int igh, double *A, double *eigenValsReal, double *eigenValsImag, int schurVectorFlag, double *eigMat);

extern double*matmul(double*A, int nA, int mA, double *B, int nB, int mB);

extern double*matsub(double*A, int nA, int mA, double *B, int nB, int mB);

void usage()
{
    printf("test_schurVectors.exe [-h | -n sizeOfMatrix | -s seed]\n");
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
    double *schurMat = (double *) calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i++)
        schurMat[i + i * n] = 1;
    // Now we call hqr. At the end schurMat will contain the schur vectors
    int ret = hqr(n,n,1,n,T,eigValsReal,eigValsImag,1,schurMat);
    if (ret != 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %d\n",ret);
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
        printf("%% [ ORTH ] n = %4d; checks = [ ", n );
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
        printf(" %1.10e", orthZ );

    //  check || A * Z - Z * T ||_F / || A ||_F
        double normR, normA;
        normR = 0.0e+00;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp = 0.0e+00;
                for (int k = 0; (k < n)&&(k < j+2); k++) {
                    tmp += schurMat0(i,k)*t0(k,j);
                }
                for (int k = 0; k < n; k++) {
                    tmp -= a0(i,k)*schurMat0(k,j);
                }
                normR += tmp * tmp ;
            }
        }
        normR = sqrt( normR );
        normA = 0.0e+00;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                normA += t0(i,j) * t0(i,j) ;
            }
        }
        normA = sqrt( normA );
        printf(" %1.10e", normR / normA );
        printf(" ];\n");
    free(A);
    free(T);
    free(schurMat);
    return 0;
}
