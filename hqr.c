#include<math.h>
#define a1(i,j) A[(i - 1) + (j - 1) * n]
#define eigenMatrix1(i,j) eigenMatrix[(i - 1) + (j - 1) * n]
extern int qrIteration(int n, double* A, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m);

extern int qrIterationVec(int n, double* A, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m, int low, int igh, double* eigenMatrix);

extern int formShift(int n, int low, double* A, int* ierr, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w);

extern int subDiagonalSearch(int n, int low, double* A, int en, double norm, double* s);

extern int doubleSubDiagonalSearch(int n, double* A, int en, int enm2, int l, double* s, double x,
        double y, double w, double* p, double* q, double* r, double* zz);

extern void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci); 

int hqr(int nm, int n, int low, int igh, double *A, double *eigenValsReal, double *eigenValsImag, int schurVectorFlag, double *eigenMatrix)
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
            return -1;
    }
postExceptionalShift_130:
    its = its + 1;
    itn = itn - 1;
    m = doubleSubDiagonalSearch(n, A, en, enm2, l, &s, x, y, w, &p, &q, &r, &zz);
    // double qr step
    if (schurVectorFlag) 
        qrIterationVec(n,A,en,na,l,&s,&x,&y,&p,&q,&r,&zz,m,low,igh,eigenMatrix);
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
            zz = eigenMatrix1(i,na);
            eigenMatrix1(i,na) = q * zz + p * eigenMatrix1(i,en);
            eigenMatrix1(i,en) = q * eigenMatrix1(i,en) - p * zz;
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
endOfProgram_1001:
    return en;
}
