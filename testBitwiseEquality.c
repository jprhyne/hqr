#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#define a0(i,j) A[(i) + (j) * n]
#define b0(i,j) B[(i) + (j) * n]
#define schurMatC0(i,j) schurMatC[(i) + (j) * n]
#define schurMatF0(i,j) schurMatF[(i) + (j) * n]

extern void hqr2_(int *nm, int *n, int *low, int *igh, double *h, double *wr,
        double *wi, double *z, int *ierr);

extern int hqr(int nm, int n, int low, int igh, double *A, double *eigenValsReal, double *eigenValsImag, int schurVectorFlag, double *eigMat);

void usage()
{
    printf("main_hqr_test.exe [-h | -n sizeOfMatrix | -s seed]\n");
    printf("\t-h: Print this help dialogue\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-s: Sets the seed. The following argument must be a positive integer.\n");
    printf("\t\tThe default value is 24875.");
}

int main (int argc, char **argv) 
{
    // Seeds the random number generator for repeatability
    int seed = 24875;
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
    // B will be a copy of A, and passed into hqr2.f. This is because we want to 
    // test if our port is exactly the same
    double *B = (double *) calloc(n*n,sizeof(double));
    double *eigValsRealC = (double *) malloc(n*n*sizeof(double));
    double *eigValsImagC = (double *) malloc(n*n*sizeof(double));
    double *eigValsRealF = (double *) malloc(n*n*sizeof(double));
    double *eigValsImagF = (double *) malloc(n*n*sizeof(double));
    // Generate A as a random matrix.
    for(int i = 0; i < n; i++) {
        int start = 0;
        if (i - 1 > 0)
            start = i - 1;
 	    for(int j = start; j < n; j++) {
            double val = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	        a0(i,j) = val; 
	        b0(i,j) = val; 
        }
    }
    // Create a matrix to store the Schur Vectors
    double *schurMatC = (double *) calloc(n*n,sizeof(double));
    double *schurMatF = (double *) calloc(n*n,sizeof(double));
    for (int i = 0; i < n; i++) {
        schurMatC[i + i * n] = 1;
        schurMatF[i + i * n] = 1;
    }
    // Now we call hqr. At the end schurMat will contain the schur vectors
    int ret = hqr(n,n,1,n,A,eigValsRealC,eigValsImagC,1,schurMatC);
    if (ret != 0) {
        // This means that hqr did not converge to at some index,
        // so we print it out and terminate execution as our Schur
        // vectors will not be correct
        printf("Did not converge at index: %d\n",ret);
        return 1;
    }
    // Call the (modified) hqr2.f function 
    int ione = 1;
    int ierr;
    hqr2_(&n,&n,&ione,&n,B,eigValsRealF,eigValsImagF,schurMatF,&ierr);
    // Now, we need to check if 
    // 1) A = B
    // 2) eigVals{Real,Imag}C = eigVals{Real,Imag}F
    // 3) schurMatC = schurMatF
    // These will all be done by looking at the relative error in the frobenius 
    // norm for 1 and 3, and the euclidean norm for 2
    
    // Checking 1
    double schurEq,tmp;
    schurEq = 0.0;
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++) {
            tmp = a0(i,j) - b0(i,j);
            schurEq += tmp * tmp;
        }
    schurEq = sqrt(schurEq);

    // Checking 2
    double eigRealEq,eigImagEq;
    eigRealEq = 0.0;
    eigImagEq = 0.0;
    for (int i = 0; i < n; i++) {
        tmp = eigValsRealC[i] - eigValsRealF[i];
        eigRealEq = tmp * tmp;
        tmp = eigValsImagC[i] - eigValsImagF[i];
        eigImagEq = tmp * tmp;
    }
    eigRealEq = sqrt(eigRealEq);
    eigImagEq = sqrt(eigImagEq);

    // Checking 3
    double zEq;
    zEq = 0.0;
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++) {
            tmp = schurMatC0(i,j) - schurMatF0(i,j);
            zEq += tmp * tmp;
        }
    zEq = sqrt(zEq);

    // Now we print these results to a file titled "bitEq.txt
    FILE *testingFile = fopen("bitEq.txt","a");
    fprintf(testingFile, "n=%d,schurDiff=%1.10e, eigValsRealDiff=%1.10e, eigValsImagDiff=%1.10e, zMatDiff=%1.10e\n",n,schurEq,eigRealEq,eigImagEq,zEq);
    fclose(testingFile);

    free(A);
    free(B);
    free(eigValsRealC);
    free(eigValsImagC);
    free(eigValsRealF);
    free(eigValsImagF);
    free(schurMatC);
    free(schurMatF);
    return 0;
}