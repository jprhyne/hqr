#include"schurToEigen.h"
#include<stdio.h>
/**
 * Input:
 *
 * low: lowest index value we are considering for our matrix
 * igh: highest index value we are considering for our matrix 
 * norm: The norm of the original matrix A
 * n: The number of rows/columns in the original matrix A
 * eigenValsReal: array of the real parts of eigenvalues computed by hqr
 * eigenValsImag: array of the complex parts of eigenvalues computed by hqr
 * T: Quasi-Schur form of the original matrix A computed by hqr
 * eigenMatrix: This starts as the schur matrix computed by hqr, then at the
 *      end, it will contain the eigenvectors.
 *
 */
void schurToEigen(int low, int igh, double norm, int n, double *eigenValsReal, double *eigenValsImag, double *T, double *eigenMatrix)
{
    FILE *testingFile = fopen("outputFileC.txt","w");
    int en,na,m,i,j,enm2;
    double p,q,w,r,zz,s,t,tst1,tst2,x,y,ra,sa,vr,vi;
    if (norm <= 0.)
        return;
    for (int nn = 1; nn <= n; nn++){
        /*
        for (int i = 1; i <= n; i++){
            for (int j = 1; j <= n; j++) {
                fprintf(testingFile, "%1.20e,", t1(i,j));
            }
            fprintf(testingFile, "%1.20e\n", t1(i,j));
        }
    fprintf(testingFile, "\n");
    */
        en = n + 1 - nn;
        p = eigenValsReal[en - 1];
        q = eigenValsImag[en - 1];
        na = en - 1;
        if (q > 0)
            continue;
        if (q == 0) {
            // This is the 600 label from fortran
            m = en;
            t1(en,en) = 1.0;
            if (na == 0)
                continue;
            for (int ii = 1; ii <= na; ii++) {
                i = en - ii;
                w = t1(i,i) - p;
                r = 0.0;
                for (j = m; j <= en; j++)
                    r = r + t1(i,j) * t1(j,en);
                if (eigenValsImag[i - 1] < 0) {
                    zz = w;
                    s = r;
                    continue;
                }
                m = i;
                if (eigenValsImag[i - 1] == 0.0) {
                    t = w;
                    if (t == 0.0) {
                        tst1 = norm;
                        t = tst1;
                        do {
                            t = 0.01 * t;
                            tst2 = norm + t;
                        } while (tst2 >= tst1);
                    }
                    t1(i,en) = -r / t;
                } else {
                    //solve real equations
                    x = t1(i,i+1);
                    y = t1(i+1,i);
                    q = (eigenValsReal[i - 1] - p) * (eigenValsReal[i - 1] - p) + eigenValsImag[i - 1] * eigenValsImag[i - 1];
                    t = (x * s - zz * r) / q;
                    t1(i,en) = t;
                    if (fabs(x) > fabs(zz))
                        t1(i+1,en) = (-r - w * t) / x;
                    else 
                        t1(i+1,en) = (-s - y * t) / zz;
                }
                //overflow control
                t = fabs(t1(i,en));
                if (t == 0) 
                    continue;
                tst1 = t;
                tst2 = tst1 + 1.0 / tst1;
                if (tst2 > tst1)
                    continue;
                for (j = i; j <= en; j++)
                    t1(j,en) = t1(j,en) / t;
            }//label 700 restarts this loop
            // end real vector
        } else {
            m = na;
            //last vector component chosen imaginary so that the eigenvector
            // matrix is triangular
            if (fabs(t1(en,na)) <= fabs(t1(na,en))) {
                double a = 0.0;
                double b = 0.0;
                cdiv(0.0, -t1(na,en),t1(na,na)-p,q,&a,&b);
                t1(na,na) = a;
                t1(na,en) = b;
            } else {
                t1(na,na) = q / t1(en,na);
                t1(na,en) = -(t1(en,en) - p) / t1(en,na);
            }
            t1(en,na) = 0.0;
            t1(en,en) = 1.0;
            enm2 = na - 1;
            if (enm2 == 0)
                continue;
            for(int ii = 1; ii <= enm2; ii++) { //fortran ends this loop at 795
                i = na - ii;
                w = t1(i,i) - p;
                ra = 0.0;
                sa = 0.0;
                for (j = m; j <= en; j++) {
                    ra = ra + t1(i,j) * t1(j,na);
                    sa = sa + t1(i,j) * t1(j,en);
                }
                if (eigenValsImag[i - 1] < 0) {
                    zz = w;
                    r = ra;
                    s = sa;
                    continue;
                }
                m = i;
                if (eigenValsImag[i - 1] == 0.0) {
                    double a = 0.0;
                    double b = 0.0;
                    cdiv(-ra, -sa, w, q, &a, &b);
                    t1(i,na) = a;
                    t1(i,en) = b;
                } else {
                    // solve complex equations
                    x = t1(i, i + 1);
                    y = t1(i + 1, i);
                    vr = (eigenValsReal[i - 1] - p) * (eigenValsReal[i - 1] - p) + eigenValsImag[i - 1] * eigenValsImag[i - 1] - q * q;
                    vi = (eigenValsReal[i - 1] - p) * 2.0 * q;
                    if (vr == 0.0 && vi == 0.0) {
                        tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
                        vr = tst1;
                        do {
                            vr = 0.01 * vr;
                            tst2 = tst1 + vr;
                        } while (tst2 > tst1);
                    }//784 is after this if block
                    double a = 0.0;
                    double b = 0.0;
                    cdiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi,  &a, &b);
                    t1(i,na) = a;
                    t1(i,en) = b;
                    if (fabs(x) > fabs(zz) + fabs(q)) {
                        t1(i + 1, na) = (-ra - w * t1(i,na) + q * t1(i,en)) / x;
                        t1(i + 1, en) = (-sa - w * t1(i,en) - q * t1(i,na)) / x;
                    } else {
                        a = 0.0;
                        b = 0.0;
                        cdiv(-r - y * t1(i, na), -s - y * t1(i,en), zz, q, &a, &b);
                        t1(i + 1, na) = a;
                        t1(i + 1, en) = b;
                    }
                }
                // overflow control 790
                t = fabs(t1(i, na));
                if (t < fabs(t1(i, en)))
                    t = fabs(t1(i, en));
                if ( t == 0.0 )
                    continue;
                tst1 = t;
                tst2 = tst1 + 1.0 / tst1;
                if (tst2 > tst1)
                    continue;
                for (j = i; j <= en; j++) {
                    t1(j,na) = t1(j,na) / t;
                    t1(j,en) = t1(j,en) / t;
                }
            } //label 795 restarts this loop
            // end complex vector
        }
    }
    // This section produces the vectors of isolated roots. 
    // In our testing we do nothing with this, however this is kept in the event
    // that balance or a similar function is ported and/or used
    for (i = 1; i <= n; i++) {
        if (i >= low && i <= igh) 
           continue;
        for (j = 1; j <= n; j++)
            eigenMatrix1(i,j) = t1(i,j);
    }
    for (int jj = low; jj <= n; jj++) {
        j = n + low - jj;
        m = j;
        if (m > igh)
            m = igh;
        for (i = low; i <= igh; i++) {
            zz = 0.0;
            for (int k = low; k <= m; k++) 
                zz = zz + eigenMatrix1(i,k) * t1(k,j);
            eigenMatrix1(i,j) = zz;
        }
        // Do not add anything here
    }
    return;
}
