#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x
#define b0(i,j) B[(i) + (j) * n]
#define a0(i,j) A[(i) + (j) * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define a1(i,j) A[(i - 1) + (j - 1) * n]

// This may be bad practice, but we put all malloc'd entities
// as global variables in order to make freeing the
// memory easier
double* A;
double* B;
double* wr;
double* wi;
double* eigenValsReal;
double* eigenValsImag;

/*
 * This is header is for the subroutine inside hqr.f
 * If recompiled, check the name of the function using
 * ```nm hqr.o``` and change below accordingly.
 * On a unix machine, something like
 * ```nm hqr.o | cut -d " " -f 3``` gives the correct output
 */
extern void hqr_( int* nm, int* n, int* low, int* igh,
	double* h, double* wr, double* wi, int* ierr);

/*
 * This function is going to be similar to the one above, however
 * we are instead copying one loop at a time
 */
extern void qrIteration(int n, double* h, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m);

extern int formShift(int n, int low, double* B, int* ierr, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w);

extern int subDiagonalSearch(int n, int low, double* B, int en, double norm, double* s);

extern int doubleSubDiagonalSearch(int n, double* B, int en, int enm2, int l, double* s, double x,
        double y, double w, double* p, double* q, double* r, double* zz);

void usage()
{
    printf("main_hqr_test.exe [-n sizeOfMatrix | -v | -h]\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-v: verbose flag that prints out the results of eispack hqr\n");
    printf("\t\tBy default, nothing is printed to the console\n");
    printf("\t-h: Print this help dialogue\n");
}

void freeMemory()
{
    free(A);
    free(B);
    free(wr);
    free(wi);
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

	// Seeds the random number generator for repeatability
	srand(734);
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
            }
	}

	// Allocate the memory for A to be generated. It will contain n^2 
	// elements where each element is a double precision floating point number
	A = (double *) malloc( n * n *  sizeof(double));
	B = (double *) malloc( n * n *  sizeof(double));
	// Create a vector to store the real parts of the eigenvalues
	wr = (double *) malloc( n *  sizeof(double));
	// Create a vector to store the real parts of the eigenvalues
	wi = (double *) malloc( n *  sizeof(double));

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
	/*
	 * Here, we print out A for finding the eigenvalues via MATLAB
	 * This must be done before calling hqr_ because it destroys A
	 */
        if (printFlag) {
            printf("A = [\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%3.8f,", a0(i,j));
                }
                printf("\n");
            }
            printf("]\n");
        }
	hqr_( &n, &n, &ione, &n, A, wr, wi, &ierr);
	/*
	 * Below prints out wr and wi for inspection via 
	 * MATLAB/visual inspection
	 */
        if (printFlag) {
            printf("wr = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", wr[i]);
            printf("]\n");
            printf("wi = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", wi[i]);
            printf("]\n");
        }

        // Now, we write our version. This next section is just
        // a C implementation of EISPACK's HQR.f where we move from an upper
        // Hessenberg matrix to the Schur form.
        
        eigenValsReal = (double *) malloc(n * sizeof(double));
        eigenValsImag = (double *) malloc(n * sizeof(double));

        int indexOfError = 0;
        double norm = 0;
        int k = 1;
        // These deal with our boundary conditions
        int low = 1;
        int igh = n;
        int en,m,mm,notLast,itn,its,na,enm2,l,ll,retVal,mp2;
        double x,y,z,t,w,s,r,q,p,zz,tst1,tst2;
        // Converting one section at a time
        // This section is not being used in our case until a version of
        // balance is ported
        for (i = 1; i <= n; i++) {
            for (j = k; j <= n; j++) {
                norm += fabs(b1(i,j));
            }
            k = i;
            if (i >= low && i <= igh)
                continue;
            eigenValsReal[i-1] = b1(i,i);
            eigenValsImag[i-1] = 0.0;
        }
        //initializing some variables
        en = igh;
        t = 0.0;
        itn = 30 * n;
beginEigSearch_60:
        if (en < low)
            goto endOfProgram_1001;
        its = 0;
        na = en - 1;
        enm2 = na -1;
subDiagonalSearch_70:
        l = subDiagonalSearch(n,low,B,en,norm,&s);
formShift_100:
        retVal = formShift(n,low, B, &ierr, its, itn, en, l, &s, &t, &x, &y, &w);
        // In order to emulate the behavior of the fortran code, instead 
        // of jumping to the right code inside there, we instead set a
        // return value and check what it is on exit
        // if retVal did not change from 0, we went through the entire
        // form shift section
        // if retVal is 1, then we found a single root
        // if retVal is 2, then we found a double root
        // if retVal is 3, then we did not converge and terminate with error
        switch (retVal) {
            case 0: 
                // full termination
                break;
            case 1: 
                // single root
                goto singleRoot_270;
            case 2:
                // double root
                goto doubleRoot_280;
            case 3:
                // Error termination
                goto errorThenEnd_1000;
            default:
                // This should never happen, so if it does we free memory
                // print an error message, then terminate.
                freeMemory();
                printf("Error in fortran subroutine. Check if assignment of retVal is correct\n");
                return 2;
        }
postExceptionalShift_130:
        its = its + 1;
        itn = itn - 1;
        m = doubleSubDiagonalSearch(n, B, en, enm2, l, &s, x, y, w, &p, &q, &r, &zz);
        // double qr step
        qrIteration(n,B,en,na,l, &s,&x,&y,&p,&q,&r,&zz,m);
        goto subDiagonalSearch_70;

singleRoot_270:
        eigenValsReal[en - 1] = x + t;
        eigenValsImag[en - 1] = 0;
        en = na;
        goto beginEigSearch_60;
doubleRoot_280:
        p = (y - x) / 2.0;
        q = p * p + w;
        zz = sqrt(fabs(q));
        x = x + t;
        if (q < 0)
            goto complexPair_320;
        // real pair
        if (p >= 0)
            zz = p + zz;
        else 
            zz = p - zz;
        eigenValsReal[na - 1] = x + zz;
        eigenValsReal[en - 1] = eigenValsReal[na - 1];
        if (zz != 0)
            eigenValsReal[en - 1] = x - w / zz;
        eigenValsImag[na - 1] = 0;
        eigenValsImag[en - 1] = 0;
        goto postDoubleRoot_330;
complexPair_320:
        eigenValsReal[na - 1] = x + p;
        eigenValsReal[en - 1] = x + p;
        eigenValsImag[na - 1] = zz;
        eigenValsImag[en - 1] = -zz;
postDoubleRoot_330:
        en = enm2;
        goto beginEigSearch_60;
errorThenEnd_1000:
        indexOfError = en; 
endOfProgram_1001:
        if (printFlag){
            printf("eigValReal = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", eigenValsReal[i]);
            printf("]\n");
            printf("eigValImag = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", eigenValsImag[i]);
            printf("]\n");
            printf("B = [\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%3.8f,", b0(i,j));
                }
                printf("\n");
            }
            printf("]\n");
			// TODO: Add 2 arrays wr-eigreal and wi-eigimag
        }
        freeMemory();
        return 0;
}
