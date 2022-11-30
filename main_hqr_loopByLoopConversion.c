#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define b0(i,j) B[(i) + (j) * n]
#define a0(i,j) A[(i) + (j) * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define a1(i,j) A[(i - 1) + (j - 1) * n]
#define z1(i,j) eigenMatrix[(i - 1) + (j - 1) * n]

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
double* eigenMatrix;

extern int qrIteration(int n, double* B, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m);

extern int qrIterationVec(int n, double* B, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m, int low, int igh, double* eigenMatrix);

extern int formShift(int n, int low, double* B, int* ierr, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w);

extern int subDiagonalSearch(int n, int low, double* B, int en, double norm, double* s);

extern int doubleSubDiagonalSearch(int n, double* B, int en, int enm2, int l, double* s, double x,
        double y, double w, double* p, double* q, double* r, double* zz);

extern void hqr2_(int *nm, int *n, int *low, int *igh, double *h, double *wr,
        double *wi, double *z, int *ierr);

extern void hqr_(int *nm, int *n, int *low, int *igh, double *h, double *wr,
        double *wi, int *ierr);

extern void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci); 

void usage()
{
    printf("main_hqr_test.exe [-n sizeOfMatrix | -v | -h | -t | -h | --jobv | --schur]\n");
    printf("\t-n: The following argument must be a positive integer\n");
    printf("\t\tThe default value is 20\n");
    printf("\t-v: verbose flag that prints out the results of eispack hqr\n");
    printf("\t\tBy default, nothing is printed to the console\n");
    printf("\t-t: testing flag that only prints the expected vs the\n");
    printf("\t\tactual computed eigenvalues.\n");
    printf("\t-h: Print this help dialogue\n");
    printf("\t--jobv: Flag that tells us to compute the eigen vectors\n");
    printf("\t\tOverrides --schur\n");
    printf("\t--schur: Flag that tells us to compute the Schur vectors\n");
    printf("\t\tOverridden by --jobv\n");
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
    free(eigenMatrix);
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
        int schurVectorFlag = 0;

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
            } else if ( strcmp ( *(argv + i), "--schur" ) == 0 ) {
                schurVectorFlag = 1;
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
	A = (double *) calloc( n * n, sizeof(double));
	B = (double *) calloc( n * n, sizeof(double));
	// Create a vector to store the real parts of the eigenvalues
	wr = (double *) malloc( n *  sizeof(double));
	// Create a vector to store the real parts of the eigenvalues
	wi = (double *) malloc( n *  sizeof(double));

        z = (double *) calloc( n * n, sizeof(double)); //This needs to start as the identity matrix
        eigenMatrix = (double *) calloc( n * n, sizeof(double));//This needs to start as the identity matrix
        // Start z and eigenMatrix as the identity matrix
        for (i = 0; i < n; i++) {
            z[i + i * n] = 1;
            eigenMatrix[i + i * n] = 1;
        }
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
        if (printFlag && !testFlag) {
            printf("A = [\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%3.8f,", a0(i,j));
                }
                printf("\n");
            }
            printf("]\n");
        }
        if (eigenVectorFlag || schurVectorFlag) {
            hqr2_(&n, &n, &ione, &n, A, wr, wi, z, &ierr);
        } else {
            hqr_(&n,&n,&ione,&n,A,wr,wi,&ierr);
        }
	/*
	 * Below prints out wr and wi for inspection via 
	 * MATLAB/visual inspection
	 */
        if (printFlag && !testFlag) {
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
        double x,y,z,t,w,s,r,q,p,zz,tst1,tst2,ra,sa,vr,vi;
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
        if (en < low) {
            if (eigenVectorFlag) {
                goto backSub_340;
            } else {
                goto endOfProgram_1001;
            }
        }
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
        if (eigenVectorFlag || schurVectorFlag) 
            qrIterationVec(n,B,en,na,l,&s,&x,&y,&p,&q,&r,&zz,m,low,igh,eigenMatrix);
        else 
            qrIteration(n,B,en,na,l, &s,&x,&y,&p,&q,&r,&zz,m);
        // For debugging purposes, we print out the contents of b1 to a file
		/*
        for (int i = 1; i <= n; i++){
            for (int j = 1; j < n; j++) {
                fprintf(testingFile, "%1.20f,", b1(i,j));
            }
            fprintf(testingFile, "%1.20f\n", b1(i,j));
        }
        fprintf(testingFile, "\n");
		*/
        goto subDiagonalSearch_70;

singleRoot_270:
        if (eigenVectorFlag) {
            b1(en,en) = x + t;
        }
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
        if (eigenVectorFlag) {
            x = b1(en,na);
            s = fabs(x) + fabs(zz);
            p = x / s;
            q = zz / s;
            r = sqrt(p*p + q*q);
            p = p/r;
            q = q/r;
//c     .......... row modification ..........
            for (j = na; j <= n; j++) {
                zz = b1(na,j);
                b1(na,j) = q * zz + p * b1(en,j);
                b1(en,j) = q * b1(en,j) - p * zz;
            }
//c     .......... column modification ..........
            for (i = 1; i <= en; i++) {
                zz = b1(i,na);
                b1(i,na) = q * zz + p * b1(i,en);
                b1(i,en) = q * b1(i,en) - p * zz;
            }
//c     .......... accumulate transformations ..........
            for (i = low; i <= igh; i++) {
                zz = z1(i,na);
                z1(i,na) = q * zz + p * z1(i,en);
                z1(i,en) = q * z1(i,en) - p * zz;
            }
        }
        goto postDoubleRoot_330;
complexPair_320:
        eigenValsReal[na - 1] = x + p;
        eigenValsReal[en - 1] = x + p;
        eigenValsImag[na - 1] = zz;
        eigenValsImag[en - 1] = -zz;
postDoubleRoot_330:
        en = enm2;
        goto beginEigSearch_60;
backSub_340:
        if (norm == 0.)
            goto endOfProgram_1001;
        for (int nn = 1; nn <= n; nn++) {
            en = n + 1 - nn;
            p = eigenValsReal[en - 1];
            q = eigenValsImag[en - 1];
            na = en -1;
                //fortran says if (q) 710, 600, 800
            if (q < 0) {
                goto complexVector_710;
            } else if (q > 0) {
                goto b800; 
            }
a600:       m = en;
            b1(en,en) = 1.0;
            if (n == 0.0) {
                goto b800;
            }
            for (int ii = 1; ii <= na; ii++) {
                i = en - ii;
                w = b1(i,i) - p;
                r = 0.0;
                for (j = m; j<= en; j++)
                    r = r + b1(i,j) * b1(j,en);
                if (eigenValsImag[i - 1] < 0.0) {
                    zz = w;
                    s = r;
                    continue;
                }
                m = i;
                if (eigenValsImag[i - 1] == 0) {
                    t = w;
                    if ( t == 0.0 ) {
                        do {
                            tst1 = norm;
                            t = tst1;
                            t = 0.01 * t;
                            tst2 = norm + t;
                        } while (tst2 > tst1);
                        b1(i,en) = -r / t;
                    }
                    goto overflowControl_680;
                }
//c			.......... solve real equations ..........
                x = b1(i,i+1);
                y = b1(i+1,i);
                q = (eigenValsReal[i - 1] - p) * (eigenValsReal[i - 1] - p) + eigenValsImag[i - 1] * eigenValsImag[i - 1];
                t = (x * s - zz * r) / q;
                b1(i,en) = t;
                if (fabs(x) > fabs(zz)) {
                    b1(i+1, en) = (-r - w * t) / x;
                } else {
                    b1(i + 1, en) = (-s - y * t) / zz;
                }
//c				.......... overflow control ..........
overflowControl_680:
                t = fabs(b1(i,en));
                if ( t == 0.0) continue; //go to 700 (end of for loop)
                tst1 = t;
                tst2 = tst1 + 1.0/tst1;
                if ( tst2 > tst1) continue; //go to 700 (end of loop)
                for (j = i; j <= en; j++)
                    b1(j,en) = b1(j,en) / t;
            }
        }
//c		.......... end real vector ..........
        goto b800;
//c		.......... complex vector ..........
complexVector_710:
        m = na;
//c     .......... last vector component chosen imaginary so that
//c                eigenvector matrix is triangular ..........
        if (fabs(b1(en,na)) <= fabs(b1(na,en))) { //go to 720
            double a;
            double b;
            cdiv(0.0,-b1(na,en),b1(na,na)-p,q,&a,&b);
            b1(na,na) = a;
            b1(na,en) = b;
        } else {
            b1(na,na) = q / b1(en,na);
            b1(na,en) = -(b1(en,en) - p) / b1(en,na);
        }
        b1(en,na) = 0.0;
        b1(en,en) = 1.0;
        enm2 = na - 1;
        if (enm2 == 0) goto b800; // go to 800
        for (int ii = 1; ii<=enm2;ii++){
            i = na - ii;
            w = b1(i,i) - p;
            ra = 0.0;
            sa = 0.0;
            for (j = m; j<=en; j++) {
                ra = ra + b1(i,j) * b1(j,na);
                sa = sa + b1(i,j) * b1(j,en);
            }
            if (eigenValsImag[i - 1] < 0.0) {
                zz = w;
                r = ra;
                s = sa;
                continue;
            }
            m = i;
            if (eigenValsImag[i - 1] == 0) {
                double a;
                double b;
                cdiv(-ra,-sa,w,q,&a,&b);
                b1(i,na) = a;
                b1(i,en) = b;
                goto overflowControl_790;
            }
//c			.......... solve complex equations 
            x = b1(i,i+1);
            y = b1(i+1,i);
            vr = (eigenValsReal[i - 1] - p) * (eigenValsReal[i-1] - p) + eigenValsImag[i-1] * eigenValsImag[i-1] - q * q;
            vi = (eigenValsReal[i-1] - p) * 2.0 * q;
            if (vr == 0.0 && vi== 0.0) {
                tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
                do {
                    vr = tst1;
                    vr = 0.01 * vr;
                    tst2 = tst1 + vr;
                } while (tst2 > tst1);
            }
            double a;
            double b;
            cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,&a,&b);
            b1(i+1,na) = a;
            b1(i+1,en) = b;
            if (fabs(x) > fabs(zz) + fabs(q)) {
                b1(i+1,na) = (-ra - w * b1(i,na) + q * b1(i,en)) / x;
                b1(i+1,en) = (-sa - w * b1(i,en) - q * b1(i,na)) / x;
            } else {
                double a;
                double b;
                cdiv(-r-y*b1(i,na),-s-y*b1(i,en),zz,q,&a,&b);
                b1(i+1,na) = a;
                b1(i+1,en) = b;
            }
overflowControl_790:
            if (fabs(b1(i,na)) >= fabs(b1(i,en))) {
                t = fabs(b1(i,na));
            } else {
                t = fabs(b1(i,en));
            }
            if (t == 0.0) continue;
            tst1 = t;
            tst2 = tst1 + 1.0/tst1;
            if (tst2 > tst1) continue;
            for (j = i; j <= en; j++) {
                b1(j,na) = b1(j,na) / t;
                b1(j,en) = b1(j,en) / t;
            }
        }	
//800 end complex vector
b800:
//c		.......... end back substitution.
//c				   vectors of isolated roots ..........
        for (i = 1; i <= n; i++) {
            if ( i >= low && i <= igh ) continue;
            for (j = i; j<=n; j++)
                z1(i,j) = b1(i,j);
        }
//c     .......... multiply by transformation matrix to give
//c                vectors of original full matrix.
//c                for j=n step -1 until low do -- ..........
        for (int jj = low; jj <= n; jj++) {
            j = n + low - jj;
            m = j;
            if (igh < j)
                m = igh;
            for (i = low; i <= igh; i++) {
                zz = 0.0;
                for (int k = low; k <= m; k++) 
                    zz = zz + z1(i,k) * b1(k,j);
                z1(i,j) = zz;
            }
        }
        goto endOfProgram_1001;
errorThenEnd_1000:
        indexOfError = en; 
endOfProgram_1001:
	for (int i = 0; i < n; i++) {
            eigRealDiff[i] = eigenValsReal[i] - wr[i];
            eigImagDiff[i] = eigenValsImag[i] - wi[i];
	}
        // Compute the 1-norms of eigRealDiff and eigImagDiff
        double normReal = 0.;
        double normImag = 0.;
        for (int i = 0; i < n; i ++) {
            if (eigRealDiff[i] != 0.)
                normReal += fabs(eigRealDiff[i]);
            if (eigImagDiff[i] != 0.)
                normImag += fabs(eigImagDiff[i]);
        }
        if (printFlag && !testFlag){
            printf("eigValReal = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", eigenValsReal[i]);
            printf("]\n");
            printf("eigValImag = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %5.8f,\n", eigenValsImag[i]);
            printf("]\n");
            printf("eigRealDiff = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %1.20f,\n", eigRealDiff[i]);
            printf("]\n");
            printf("eigImagDiff = [\n");
            for ( int i = 0; i < n; i++)
                printf("    %1.20f,\n", eigImagDiff[i]);
            printf("]\n");
            printf("B = [\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%3.8f,", b0(i,j));
                }
                printf("\n");
            }
            printf("]\n");
        } else if (testFlag) {
            printf("Seed=%d, diff=%1.20f, diff is 0: %d\n",seed,normReal + normImag,normReal + normImag == 0.0);
		}
        freeMemory();
        return 0;
}
