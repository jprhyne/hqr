#include "hqr.h"
/**
 * nm: leading dimension for A
 * n: number of columns in A
 * low: lowest index value we are considering for our matrix
 * igh: highest index value we are considering for our matrix 
 * A: The upper hessenberg matrix we are finding the eigenvalues of
 *      and, and computing the Quasi-Schur form of
 * eigenValsReal: Array to store the real parts of the eigenvalues of A
 * eigenValsImag: Array to store the imaginary parts of the eigenvalues of A
 * schurVectorFlag: flag that determines of we compute the schur vectors or not.
 *      if we are NOT computing the schur vectors, then schurMatrix can be a null ptr
 * schurMatrix: 2D array that will store the SchurVectors. This must start as either
 *      the Identity matrix or as a collection of transformations such that
 *      B = schurMatrix * A * schurMatrix^T (ie a hessenberg reduction to A)
 *
 * return:
 *  The return value takes 3 possible ranges of values: -n-1, [-n,0), and [0,\inf)
 *  -n-1: This happens only when there is an issue in the call to formShift where 
 *      a value other than 0,1,2,3 is returned. If this happens, make sure
 *      to check the return values of formShift
 *  [-n,0): then it will be an integer value containing the negative index of 
 *      where we failed to converge. IE: -2 means that we failed to converge 
 *      when we were trying to compute the 2nd eigenvalue. So the eigenvalues 
 *      3,...,n are correct
 *  [0,\inf): Then it will be a floating point number containing the sum of
 *      the absolute values of the upper triangular part of A
 */
double hqr(int nm, int n, int low, int igh, double *A, double *eigenValsReal, double *eigenValsImag, int schurVectorFlag, double *schurMatrix)
{

    int indexOfError = 0;
    double norm = 0.0;
    // These deal with our boundary conditions
    int i,j,en,m,mm,notLast,itn,its,l,ll,retVal,mp2;
    double x,y,z,t,w,s,r,q,p,zz,tst1,tst2,ra,sa,vr,vi;
    // This section looks for any isolated eigenvalues. This is
    // really only useful if a function like 'balance' is ported or used 
    for (i = 0; i < n; i++) {
        for (j = (i == 0) ? (i):(i-1); j <= n; j++) {
            norm += fabs(a0(i,j));
        }
        if (i >= low && i <= igh)
            continue;
        eigenValsReal[i] = a0(i,i);
        eigenValsImag[i] = 0.0;
    }
    //initializing some variables
    en = igh;
    t = 0.0;
    itn = 30 * n;
    // This flag tells us if we did not do a QR Step, and only if we did, do we 
    // reset the iterations. This means we found an eigenvalue
    // and thus need to start looking for the next one
    int didQRStep = 0;
    while (en >= low) {
        if (!didQRStep) 
            its = 0;
        didQRStep = 0; //Resetting after we use its value to prevent infinite loops
        l = subDiagonalSearch(n,low,A,en,norm,&s);
        // Formshift does not need updating, just to change the inputs
        retVal = formShift(n,low, A, its, itn, en, l, &s, &t, &x, &y, &w);
        // In order to emulate the behavior of the fortran code, instead 
        // of jumping to the right code inside there, we instead set a
        // return value and check what it is on exit
        // if retVal did not change from 0, we went through the entire
        // form shift section
        // if retVal is 1, then we found a single root
        // if retVal is 2, then we found a double root
        // if retVal is 3, then we did not converge and terminate with error
        // Any other value means there was an error inside formShift, and is
        // not supported by this implementation
        switch (retVal) {
            case 0: 
                // full termination
                // This means we need to perform our QR step then look again for an eigenvalue
                its = its + 1;
                itn = itn - 1;
                m = doubleSubDiagonalSearch(n, A, en, l, &s, x, y, w, &p, &q, &r, &zz);
                // double qr step
                qrIteration(n,A,en,l,&s,&x,&y,&p,&q,&r,&zz,m,schurVectorFlag,low,igh,schurMatrix);
                didQRStep = 1;
                break;
            case 1: 
                // single root was found
                if (schurVectorFlag)
                    a0(en,en) = x + t;
                eigenValsReal[en] = x + t;
                eigenValsImag[en] = 0.0;
                en = en - 1;
                break;
            case 2:
                // double root
                p = (y - x) / 2.0;
                q = p * p + w;
                zz = sqrt(fabs(q));
                if (schurVectorFlag) {
                    a0(en,en) = x + t;
                    a0(en - 1,en - 1) = y + t;
                }
                x = x + t;
                if (q < 0) {
                    // Complex pair
                    eigenValsReal[en - 1] = x + p;
                    eigenValsReal[en] = x + p;
                    eigenValsImag[en - 1] = zz;
                    eigenValsImag[en] = -zz;
                } else {
                    // real pair
                    if (p >= 0)
                        zz = p + zz;
                    else 
                        zz = p - zz;
                    eigenValsReal[en - 1] = x + zz;
                    eigenValsReal[en] = eigenValsReal[en - 1];
                    if (zz != 0)
                        eigenValsReal[en] = x - w / zz;
                    eigenValsImag[en - 1] = 0.0;
                    eigenValsImag[en] = 0.0;

                    if (schurVectorFlag) {
                        x = a0(en,en - 1);
                        s = fabs(x) + fabs(zz);
                        p = x / s;
                        q = zz / s;
                        r = sqrt(p*p + q*q);
                        p = p/r;
                        q = q/r;
                //c     .......... row modification ..........
                        for (j = en - 1; j < n; j++) {
                            zz = a0(en - 1, j);
                            a0(en - 1, j) = q * zz + p * a0(en, j);
                            a0(en, j) = q * a0(en, j) - p * zz;
                        }
                //c     .......... column modification ..........
                        for (i = 0; i <= en; i++) {
                            zz = a0(i, en - 1);
                            a0(i, en - 1) = q * zz + p * a0(i, en);
                            a0(i, en) = q * a0(i, en) - p * zz;
                        }
                //c     .......... accumulate transformations ..........
                        for (i = low; i <= igh; i++) {
                            zz = schurMatrix0(i, en - 1);
                            schurMatrix0(i, en - 1) = q * zz + p * schurMatrix0(i, en);
                            schurMatrix0(i, en) = q * schurMatrix0(i, en) - p * zz;
                        }
                    }
                }
                en = en - 2;
                break;
            case 3:
                // Error termination, we now return the index at which we failed to find the eigenvalue
                return -en;
            default:
                // This should never happen, however if it does through an error
                // in formShift, we return a nonsense value.
                return -n - 1;
        }
    }
    return norm;
}
