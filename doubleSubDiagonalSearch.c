#define b0(i,j) B[(i) + (j) * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#include<math.h>
/**
 *  This file is responsible for the form shift behavior. 
 *  The intention is that it will be a translation of 
 *  after the modification of itn,its from around label 130  
 *  to label 160 
 */
/**
 *  Inputs that are passed as a pointer are modified after calling
 *  this function. 
 *  the index m is returned after this 
 */
int doubleSubDiagonalSearch(int n, double* B, int en, int enm2, int l, double* s, 
        double x, double y, double w, double* p, double* q, double* r, double* zz,
        int eigenVectorFlag)
{
    int m;
    for (int mm = l; mm <= enm2; mm++) {
        m = enm2 + l - mm;
        *zz = b1(m,m);
        *r = x - *zz;
        *s = y - *zz;
        *p = (*r * *s - w) / b1(m+1,m) + b1(m,m+1);
        *q = b1(m+1,m+1) - *zz - *r - *s;
        *r = b1(m+2,m+1);
        *s = fabs(*p) + fabs(*q) + fabs(*r);
        *p = *p / *s;
        *q = *q / *s;
        *r = *r / *s;
        if (m == l) break;
        double tst1 = fabs(*p) * (fabs(b1(m-1, m-1)) + fabs(*zz) + fabs(b1(m+1,m+1)));
        double tst2 = tst1 + fabs(b1(m,m-1)) * (fabs(*q) + fabs(*r));
        if (tst1 == tst2) break;
    }
    for (int i = m + 2; i <= en; i++) {
        b1(i,i-2) = 0.0;
        if (i == m + 2) continue;
        b1(i,i-3) = 0.0;
    }
    return m;
}
