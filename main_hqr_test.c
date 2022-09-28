#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x

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
    printf("main_hqr_test.exe [-n sizeOfMatrix, -v]\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-v: verbose flag that prints out the results of eispack hqr\n");
    printf("\t\tBy default, nothing is printed to the console\n");
    printf("\t-h: Print this help dialogue\n");
}

void freeMemory(void *A, void *B, void*wr, void*wi){
    free(A);
    free(B);
    free(wr);
    free(wi);
}

int main(int argc, char ** argv) {
	/*
	 * Declaring variables to pass into the hqr subroutine
	 */
	int i, ierr, j, n;
	double *A, *B, *wr, *wi;
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
			A[i+j*n] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	
	// since hqr needs A to be upper hessenberg, we will set all other values 
	// to 0	in order find the eigenvalues using other sources
	for ( int i = 2; i < n; i++ ) {
		for (int j = 0; j < i - 1; j++) {
			A[i + j*n] = 0;
		}
	}
        // Store a copy of A into B so that we can run our own
        // version of hqr written in C 
        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                B[i +j*n] = A[i + j*n];
	/*
	 * Here, we print out A for finding the eigenvalues via MATLAB
	 * This must be done before calling hqr_ because it destroys A
	 */
        if (printFlag) {
            printf("A = [\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n - 1; j++) {
                    printf("%f,", A[i+j*n]);
                }
                printf("%f\n", A[i + (n-1)*n]);
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
        
        double *eigenValsReal, *eigenValsImag;
        eigenValsReal = (double *) malloc(n * sizeof(double));
        eigenValsImag = (double *) malloc(n * sizeof(double));

        int indexOfError = 0;
        double norm = 0;
        int k = 0;

        int low = 0;
        int igh = n-1;

        for (int i = 0; i < n; i++ ){
            for (int j = k; j < n; j++) 
                norm += abs(B[i + j * n]); // 1 norm of the matrix A
            k = i;
            // This never executes, however it is functionality in the 
            // original hwr subroutine, so I am repeating it here
            // in case there is a chance of balanc being reimplemented as well.
            if ( i <= low || i >= igh) {
                wr[i] = B[i + i*n];
                wi[i] = 0;
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

            int l = 0;
            // Look for single small sub-diagonal element 
            for (int i = low; i <= en; i++) {
                l = en + low - i;
                if (l != low) {
                    double s = abs(B[(l-1) + (l-1)*n]) + abs(B[l + l*n]);
                    if (s == 0) //consider using another method of checking if s is 0 
                        s = norm;
                    double tst1 = s;
                    double tst2 = tst1 + abs(abs(B[l + (l-1)*n]));
                    // Found it?
                    if (tst1 == tst2)
                        break;
                }
            }
            // Form shift
            x = B[en + en*n];
            // Single root found. Why?
            if (l == en) {
                wr[en] = x + t;
                wi[en] = 0;
                en = na;
                continue;
            }
            y = B[na + na*n];
            w = B[en + na*n] * B[na + en*n];

            //Double root found. Why?
            if (l == na) {
                p = (y-x) / 2;
                q = p*p+w;
                zz = sqrt(abs(q));
                x += t;
                if (q >= 0) {
                    // Real pair
                    zz = (p >= 0)?(p + zz):(p - zz);
                    wr[na] = x + zz;
                    wr[en] = (zz == 0)?(wr[na]):(x - w / zz);
                    wi[na] = 0;
                    wi[en] = 0;
                } else {
                    // Complex pair
                    wr[na] = x + p;
                    wr[en] = x + p;
                    wi[na] = zz;
                    wi[en] = -zz;
                }
            } else if (itn == 0){
                // We have not converged to all eigenvalues after 30*n iterations.
                ierr = en;
                // To prevent memory leaks
                freeMemory(A,B,wr,wi);
                return 1;
            } else if (its == 10 || its == 20) {
                // form exceptional shift
                t += x;
                for (int j = low; j <= en; j++) 
                    B[i+i*n] -= x;
                s = abs(B[en + na*n]) + abs(B[na + enm2*n]);
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
                zz = B[m+m*n];
                r = x - zz;
                s = y - zz;
                p = (r*s-w) / B[(m+1) + m*n] + B[m + (m+1)*n];
                q = B[m+1 + (m+1)*n] - zz - r - s;
                r = B[m+2 + (m+1)*n];
                s = abs(p) + abs(q) + abs(r);
                p /= s;
                q /= s;
                r /= s;
                if (m != l) {
                    tst1 = abs(p)*(abs(B[(m-1) + (m-1)*n]) + abs(zz) + abs(B[m+1 + (m+1)*n]));
                    tst2 = tst1 + abs(B[m + (m-1)*n])*(abs(q) + abs(r));
                    // Why not just test if right part of tst2 is 0?
                    if (tst1 == tst2)
                        break;
                }
            }
            for (int i = m-2;i <= en;i++) {
                B[i + (i-2)*n] = 0;
                if (i != m - 2)
                    B[i + (i-3)*n] = 0;
            }
            // double qr step involving rows l to en and columns m to en
            for (int k = m; k <= na; k++) {
                int lastIter = (na - k == 0)?(0):(1);
                if ( k != m ) {
                    p = B[k + (k-1)*n];
                    q = B[k+1 + (k-1)*n];
                    r = (!lastIter) ? (B[(k+2) + (k-1)*n]):(0);
                    x = abs(p) + abs(q) + abs(r);
                    if ( x == 0 )
                        continue;
                    p /= x;
                    q /= x;
                    r /= x;
                    s = sqrt(p*p + q*q + r*r);
                }
                if (p < 0)
                    s *= -1;
                if ( k != m ) {
                    B[k + (k-1)*n] = x * -s;
                } else if (l != m){
                    B[k + (k-1)*n] = -B[k + (k-1)*n];
                }
                p += s;
                x = p/s;
                y = q/s;
                zz=r/s;
                q /= p;
                r /= p;
                // Row modification
                for (int j=k; j<=en; j++) {
                        p = B[k+j*n] + q * B[k+1+j*n];
                        B[k+j*n] -= p*x;
                        B[k+1+j*n] -= p*y;
                        if (!lastIter)
                            B[k+2+j*n] -= p*y*zz;
                }
                    // Column modification
                int colModMax = (en < k+3) ? (en):(k+3);
                for (int i = l; i<= colModMax;i++) {
                    p = x*B[i+k*n] + y * B[i + (k+1)*n];
                    B[i+k*n] -= p;
                    B[i+(k+1)*n] -= p*q;
                    if (!lastIter)
                        B[i+(k+2)*n] -= p*r;
                }

            }
        }

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

        //Here, we are completely done with computing the eigenvalues.
        //Now, we free the memory and then exit

        // Freeing memory
        freeMemory(A,B,wr,wi);
	return 0;
}
