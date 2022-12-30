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
 *      if we are not computing the schur vectors, then schurMatrix can be a null ptr
 * schurMatrix: 2D array that will store the SchurVectors. This must start as either
 *      the Identity matrix or as a collection of transformations such that
 *      B = schurMatrix * A * schurMatrix^T (ie a hessenberg reduction to A)
 *
 * return:
 *  The return value takes 2 possible ranges of values: [-n,0) and [0,\inf)
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
    double norm = 0;
    int k = 1;
    // These deal with our boundary conditions
    int i,j,ierr,en,m,mm,notLast,itn,its,na,enm2,l,ll,retVal,mp2;
    double x,y,z,t,w,s,r,q,p,zz,tst1,tst2,ra,sa,vr,vi;
    // Converting one section at a time
    // This section is not being used in our case until a version of
    // balance is ported
    for (i = 1; i <= n; i++) {
        for (j = k; j <= n; j++) {
            norm += fabs(a1(i,j));
        }
        k = i;
        if (i >= low && i <= igh)
            continue;
        eigenValsReal[i-1] = a1(i,i);
        eigenValsImag[i-1] = 0.0;
    }
    //initializing some variables
    en = igh;
    t = 0.0;
    itn = 30 * n;
beginEigSearch_60:
    if (en < low) {
        goto endOfProgram_1001; // Otherwise, we are done
    }
    its = 0;
    na = en - 1;
    enm2 = na - 1;
subDiagonalSearch_70:
    l = subDiagonalSearch(n,low,A,en,norm,&s);
formShift_100:
    retVal = formShift(n,low, A, &ierr, its, itn, en, l, &s, &t, &x, &y, &w);
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
            return -n - 1;
    }
postExceptionalShift_130:
    its = its + 1;
    itn = itn - 1;
    m = doubleSubDiagonalSearch(n, A, en, enm2, l, &s, x, y, w, &p, &q, &r, &zz);
    // double qr step
    if (schurVectorFlag) 
        qrIterationVec(n,A,en,na,l,&s,&x,&y,&p,&q,&r,&zz,m,low,igh,schurMatrix);
    else 
        qrIteration(n,A,en,na,l, &s,&x,&y,&p,&q,&r,&zz,m);
    goto subDiagonalSearch_70;
singleRoot_270:
    if (schurVectorFlag)
        a1(en,en) = x + t;
    eigenValsReal[en - 1] = x + t;
    eigenValsImag[en - 1] = 0;
    en = na;
    goto beginEigSearch_60;
doubleRoot_280:
    p = (y - x) / 2.0;
    q = p * p + w;
    zz = sqrt(fabs(q));
    if (schurVectorFlag) {
        a1(en,en) = x + t;
        a1(na,na) = y + t;
    }
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

    if (schurVectorFlag) {
        x = a1(en,na);
        s = fabs(x) + fabs(zz);
        p = x / s;
        q = zz / s;
        r = sqrt(p*p + q*q);
        p = p/r;
        q = q/r;
//c     .......... row modification ..........
        for (j = na; j <= n; j++) {
            zz = a1(na,j);
            a1(na,j) = q * zz + p * a1(en,j);
            a1(en,j) = q * a1(en,j) - p * zz;
        }
//c     .......... column modification ..........
        for (i = 1; i <= en; i++) {
            zz = a1(i,na);
            a1(i,na) = q * zz + p * a1(i,en);
            a1(i,en) = q * a1(i,en) - p * zz;
        }
//c     .......... accumulate transformations ..........
        for (i = low; i <= igh; i++) {
            zz = schurMatrix1(i,na);
            schurMatrix1(i,na) = q * zz + p * schurMatrix1(i,en);
            schurMatrix1(i,en) = q * schurMatrix1(i,en) - p * zz;
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
errorThenEnd_1000:
    return -en;
endOfProgram_1001:
    return norm;
}
