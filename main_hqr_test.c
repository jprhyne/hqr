#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x
#define b0(i,j) B[i + j * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define a0(i,j) A[i + j * n]
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

void usage()
{
    printf("main_hqr_test.exe [-n sizeOfMatrix | -v | -h]\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-v: verbose flag that prints out the results of eispack hqr\n");
    printf("\t\tBy default, nothing is printed to the console\n");
    printf("\t-h: Print this help dialogue\n");
}

void freeMemory(){
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
	srand(0);
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
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			a0(i,j) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	
	// since hqr needs A to be upper hessenberg, we will set all other values 
	// to 0	in order find the eigenvalues using other sources
	for ( int i = 2; i < n; i++ ) {
		for (int j = 0; j < i - 1; j++) {
                    a0(i,j) = 0;
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
                for (int j = 0; j < n - 1; j++) {
                    printf("%f,", a0(i,j));
                }
                printf("%f\n", a0(i,n-1));
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

        int indexOfError = 1;
        double norm = 0;
        int k = 1;

        // These deal with our boundary conditions
        int low = 1;
        int igh = n;

        for (int i = 1; i <= n; i++ ){
            for (int j = k; j <= n; j++) 
                norm += abs(b1(i,j)); // 1 norm of the matrix A
            k = i;
            // This never executes, however it is functionality in the 
            // original hqr subroutine, so I am repeating it here
            // in case there is a chance of balanc being reimplemented as well.
            if ( i <= low || i >= igh) {
                eigenValsReal[i] = b1(i,i);
                eigenValsImag[i] = 0;
            }
        }
        int en = igh;
        double t = 0;
        int itn = 30*n;
        double s,x,y,w,zz,r,p,q,tst1,tst2;
        /*
         * This loop searches for the remaining (all in our case)
         * eigenvalues
         */
        while (en >= low) {
            int its = 0;
            //en and na are indices to add the current eigenValue to it's list
            int na = en -1;
            int enm2 = na -1;

            // Another index counter
            int l = 1;
            // Look for single small sub-diagonal element 
            for (int ll = low; ll <= en; ll++) {
                l = en + low - ll;
                if (l != low) {
                    double s = abs(b1(l-1,l-1)) + abs(b1(l,l));
                    if (s == 0) //consider using another method of checking if s is 0 
                        s = norm;
                    double tst1 = s;
                    double tst2 = tst1 + abs(abs(b1(l,l-1)));
                    // Found it?
                    if (tst1 == tst2)
                        break;
                }
            }
            // Form shift
            x = b1(en,en);
            // Single root found. Why?
            if (l == en) {
                eigenValsReal[en] = x + t;
                eigenValsImag[en] = 0;
                en = na;
                continue;
            }
            y = b1(na,na);
            w = b1(en,na) * b1(na,en);

            //Double root found. Why?
            if (l == na) {
                p = (y-x) / 2;
                q = p*p+w;
                zz = sqrt(abs(q));
                x += t;
                if (q >= 0) {
                    // Real pair
                    zz = (p >= 0)?(p + zz):(p - zz);
                    eigenValsReal[na] = x + zz;
                    eigenValsReal[en] = (zz == 0)?(eigenValsReal[na]):(x - w / zz);
                    eigenValsImag[na] = 0;
                    eigenValsImag[en] = 0;
                } else {
                    // Complex pair
                    eigenValsReal[na] = x + p;
                    eigenValsReal[en] = x + p;
                    eigenValsImag[na] = zz;
                    eigenValsImag[en] = -zz;
                }
            } else if (itn == 0){
                // We have not converged to all eigenvalues after 30*n iterations.
                ierr = en;
                // To prevent memory leaks
                freeMemory();
                return 1;
            } else if (its == 10 || its == 20) {
                // form exceptional shift
                t += x;
                for (int j = low; j <= en; j++) 
                    b1(i,i) -= x;
                s = abs(b1(en,na)) + abs(b1(na,enm2));
                x = s * 0.75;
                y = x;
                w = -.4375 * s * s;
            }
            its++;
            itn--;
            // Look for two consecutive small sub-diagonal elements.
            // for m=en-2 to l in step of -1
            int m;
            for (m = enm2; i >= l; i--) {
                zz = b1(m,m);
                r = x - zz;
                s = y - zz;
                p = (r * s - w) / b1(m + 1, m) + b1(m, m + 1);
                q = b1(m + 1, m + 1) - zz - r - s;
                r = b1(m + 2, m + 1);
                s = abs(p) + abs(q) + abs(r);
                p /= s;
                q /= s;
                r /= s;
                // We only do the following checks if we are not on the 
                // last iteration of this loop
                if (m != l) { 
                    tst1 = abs(p)*(abs(b1(m - 1, m - 1)) + abs(zz) + abs(b1(m + 1, m + 1)));
                    tst2 = tst1 + abs(b1(m, m - 1))*(abs(q) + abs(r));
                    // Why not just test if right part of tst2 is 0?
                    if (tst1 == tst2)
                        break;
                } 
            }
            for (int i = m-2;i <= en;i++) {
                b1(i, i - 2) = 0;
                if (i != m - 2)
                    b1(i, i - 3) = 0;
            }
            // double qr step involving rows l to en and columns m to en
            for (int k = m; k <= na; k++) {
                // If na == k, then we are on our
                // last iteration, and use this to determine
                // if we do some actions or not
                int lastIter = (na - k == 0)?(0):(1);
                // As long as we are not on our first iteration
                if ( k != m ) {
                    p = b1(k, k - 1);
                    q = b1(k + 1, k - 1);
                    r = (!lastIter) ? (b1( k + 2, k - 1)):(0);
                    x = abs(p) + abs(q) + abs(r);
                    if ( x == 0 )
                        continue;
                    p /= x;
                    q /= x;
                    r /= x;
                }
                s = (p < 0) ? (-sqrt(p*p + q*q + r*r)) : (sqrt(p*p + q*q + r*r));
                if ( k != m ) {
                    b1(k, k - 1) = x * -s;
                } else if (l != m){
                    b1(k, k - 1) *= -1;
                }
                p += s;
                x = p/s;
                y = q/s;
                zz=r/s;
                q /= p;
                r /= p;
                // Row modification
                for (int j=k; j<=en; j++) {
                    if (lastIter)
                        p = b1(k, j) + q * b1(k + 1, j);
                    else
                        p = b1(k, j) + q * b1(k + 1, j) + r * b1(k + 2, j);
                    b1(k,j) -= p*x;
                    b1(k + 1, j) -= p*y;
                    if (!lastIter)
                        b1(k + 1, j) -= p * zz;
                }
                // Column modification
                int colModMax = (en < k+3) ? (en):(k+3);
                for (int i = l; i<= colModMax;i++) {
                    if (lastIter)
                        p = x * b1(i, k) + y * b1(i, k + 1);
                    else
                        p = x * b1(i, k) + y * b1(i, k + 1) + zz * b1(i , k + 2);
                    b1(i, k) -= p;
                    b1(i, k + 1) -= p*q;
                    if (!lastIter)
                        b1(i, k + 2) -= p*r;
                }

            }
        }

        if (printFlag) {
            printf("eigValsReal = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", eigenValsReal[i]);
            printf("]\n");
            printf("eigValsImag = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", eigenValsImag[i]);
            printf("]\n");
        }

        //Here, we are completely done with computing the eigenvalues.
        //Now, we free the memory and then exit

        // Freeing memory
        freeMemory();
	return 0;
}
