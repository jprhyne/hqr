#include"schurToEigen.h"
#include<stdio.h>
/**
 * Input:
 *
 * low: lowest index value we are considering for our matrix
 * igh: highest index value we are considering for our matrix 
 * norm: The norm of the origien - 1l matrix A
 * n: The number of rows/columns in the origien - 1l matrix A
 * eigenValsReal: array of the real parts of eigenvalues computed by hqr
 * eigenValsImag: array of the complex parts of eigenvalues computed by hqr
 * T: Quasi-Schur form of the origien - 1l matrix A computed by hqr
 * eigenMatrix: This starts as the schur matrix computed by hqr, then at the
 *      end, it will contain the eigenvectors.
 *
 */
void schurToEigen(int low, int igh, double norm, int n, double *eigenValsReal, double *eigenValsImag, double *T, double *eigenMatrix)
{
    int en,m;
    double p,q,w,r,zz,s,t,tst0,tst2,x,y,ra,sa,vr,vi;
    if (norm <= 0.)
        return;
    for (en = n - 1; en >= 0; en--){
        p = eigenValsReal[en];
        q = eigenValsImag[en];
        if (q > 0)
            continue;
        if (q == 0) {
            // This is the 600 label from fortran
            m = en;
            t0(en,en) = 1.0;
            if (en == 0)
                continue;
            for (int i = en - 1; i >= 0; i--) {
                w = t0(i,i) - p;
                r = 0.0;
                for (int j = m; j <= en; j++)
                    r = r + t0(i,j) * t0(j,en);
                if (eigenValsImag[i] < 0) {
                    zz = w;
                    s = r;
                    continue;
                }
                m = i;
                if (eigenValsImag[i] == 0.0) {
                    t = w;
                    if (t == 0.0) {
                        tst0 = norm;
                        t = tst0;
                        do {
                            t = 0.01 * t;
                            tst2 = norm + t;
                        } while (tst2 >= tst0);
                    }
                    t0(i,en) = -r / t;
                } else {
                    //solve real equations
                    x = t0(i,i+1);
                    y = t0(i+1,i);
                    q = (eigenValsReal[i] - p) * (eigenValsReal[i] - p) + eigenValsImag[i] * eigenValsImag[i];
                    t = (x * s - zz * r) / q;
                    t0(i,en) = t;
                    if (fabs(x) > fabs(zz))
                        t0(i+1,en) = (-r - w * t) / x;
                    else 
                        t0(i+1,en) = (-s - y * t) / zz;
                }
                //overflow control
                t = fabs(t0(i,en));
                if (t == 0) 
                    continue;
                tst0 = t;
                tst2 = tst0 + 1.0 / tst0;
                if (tst2 > tst0)
                    continue;
                for (int j = i; j <= en; j++)
                    t0(j,en) = t0(j,en) / t;
            }//label 700 restarts this loop
            // end real vector
        } else {
            m = en - 1;
            //last vector component chosen imaginary so that the eigenvector
            // matrix is triangular
            if (fabs(t0(en,en - 1)) <= fabs(t0(en - 1,en))) {
                double a = 0.0;
                double b = 0.0;
                cdiv(0.0, -t0(en - 1,en),t0(en - 1,en - 1)-p,q,&a,&b);
                t0(en - 1,en - 1) = a;
                t0(en - 1,en) = b;
            } else {
                t0(en - 1,en - 1) = q / t0(en,en - 1);
                t0(en - 1,en) = -(t0(en,en) - p) / t0(en,en - 1);
            }
            t0(en,en - 1) = 0.0;
            t0(en,en) = 1.0;
            if (en == 1)
                continue;
            for(int i = en - 2; i >= 0; i--) { //fortran ends this loop at 795
                w = t0(i,i) - p;
                ra = 0.0;
                sa = 0.0;
                for (int j = m; j <= en; j++) {
                    ra = ra + t0(i,j) * t0(j,en - 1);
                    sa = sa + t0(i,j) * t0(j,en);
                }
                if (eigenValsImag[i] < 0) {
                    zz = w;
                    r = ra;
                    s = sa;
                    continue;
                }
                m = i;
                if (eigenValsImag[i] == 0.0) {
                    double a = 0.0;
                    double b = 0.0;
                    cdiv(-ra, -sa, w, q, &a, &b);
                    t0(i,en - 1) = a;
                    t0(i,en) = b;
                } else {
                    // solve complex equations
                    x = t0(i, i + 1);
                    y = t0(i + 1, i);
                    vr = (eigenValsReal[i] - p) * (eigenValsReal[i] - p) + eigenValsImag[i] * eigenValsImag[i] - q * q;
                    vi = (eigenValsReal[i] - p) * 2.0 * q;
                    if (vr == 0.0 && vi == 0.0) {
                        tst0 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
                        vr = tst0;
                        do {
                            vr = 0.01 * vr;
                            tst2 = tst0 + vr;
                        } while (tst2 > tst0);
                    }//784 is after this if block
                    double a = 0.0;
                    double b = 0.0;
                    cdiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi,  &a, &b);
                    t0(i,en - 1) = a;
                    t0(i,en) = b;
                    if (fabs(x) > fabs(zz) + fabs(q)) {
                        t0(i + 1, en - 1) = (-ra - w * t0(i,en - 1) + q * t0(i,en)) / x;
                        t0(i + 1, en) = (-sa - w * t0(i,en) - q * t0(i,en - 1)) / x;
                    } else {
                        a = 0.0;
                        b = 0.0;
                        cdiv(-r - y * t0(i, en - 1), -s - y * t0(i,en), zz, q, &a, &b);
                        t0(i + 1, en - 1) = a;
                        t0(i + 1, en) = b;
                    }
                }
                // overflow control 790
                t = fabs(t0(i, en - 1));
                if (t < fabs(t0(i, en)))
                    t = fabs(t0(i, en));
                if ( t == 0.0 )
                    continue;
                tst0 = t;
                tst2 = tst0 + 1.0 / tst0;
                if (tst2 > tst0)
                    continue;
                for (int j = i; j <= en; j++) {
                    t0(j,en - 1) = t0(j,en - 1) / t;
                    t0(j,en) = t0(j,en) / t;
                }
            } //label 795 restarts this loop
            // end complex vector
        }
    }
    // This section produces the vectors of isolated roots. 
    // In our testing we do nothing with this, however this is kept in the event
    // that balance or a similar function is ported and/or used
    for (int i = 0; i < n; i++) {
        if (i >= low && i <= igh) 
           continue;
        for (int j = i; j < n; j++)
            eigenMatrix0(i,j) = t0(i,j);
    }
    for (int j = n - 1; j >= low; j--) {
        m = j;
        if (m > igh)
            m = igh;
        for (int i = low; i <= igh; i++) {
            zz = 0.0;
            for (int k = low; k <= m; k++) 
                zz = zz + eigenMatrix0(i,k) * t0(k,j);
            eigenMatrix0(i,j) = zz;
        }
        // Do not add anything here
    }
    return;
}
