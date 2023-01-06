#define b0(i,j) B[(i) + (j) * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#include<stdbool.h>
#include<math.h>
/**
 * This file is responsible for the QR Step 
 * n: One dimension 
 */
void qrIteration(int n, double* B, int en, int na, int l, double* s, double* x, 
        double * y, double* p, double* q, double* r, double* zz, int m)
{
    int i,j,k;
    bool notLast;
    for (k = m; k <= na; k++) {
        notLast = k != na;
        if (k != m) {
            *p = b0(k,k-1);
            *q = b0(k+1,k-1);
            *r = 0.0;
            if (notLast) {
                *r = b0(k+2,k-1);
            }
            *x = fabs(*p) + fabs(*q) + fabs(*r);
            if (*x == 0.0) {
                continue;
            }
            *p = *p / *x;
            *q = *q / *x;
            *r = *r / *x;
        }
            double insideSqrt = *p * *p + *q * *q + *r * *r;
            if (*p >= 0.0) {
                *s = sqrt(insideSqrt);
            } else {
                *s = -sqrt(insideSqrt);
            }
            if (k != m) {
                b0(k,k-1) = -*s * *x;
            } else if (l != m) {
                b0(k,k-1) = -b0(k,k-1);
            }
            *p = *p + *s;
            *x = *p / *s;
            *y = *q / *s;
            *zz = *r / *s;
            *q = *q / *p;
            *r = *r / *p;
            if (notLast) {
                // Try doing nothing depending on behavior
//c     .......... row modification ..........
                for (j = k; j <= en; j++) {
                    *p = b0(k,j) + *q * b0(k+1,j) + *r * b0(k+2,j);
                    b0(k,j) = b0(k,j) - *p * *x;
                    b0(k+1,j) = b0(k+1,j) - *p * *y;
                    b0(k+2,j) = b0(k+2,j) - *p * *zz;
                }
                if (en <= k+3) {
                    j = en;
                } else {
                    j = k+3;
                }
//c     .......... column modification ..........
                for (i = l; i <= j; i++) {
                    *p = *x * b0(i,k) + *y * b0(i,k+1) + *zz * b0(i,k+2);
                    b0(i,k) = b0(i,k) - *p;
                    b0(i,k+1) = b0(i,k+1) - *p * *q;
                    b0(i,k+2) = b0(i,k+2) - *p * *r;
                }
            } else {
//c     .......... row modification ..........
                for (j = k; j <= en; j++) {
                    *p = b0(k,j) + *q * b0(k+1,j);
                    b0(k,j) = b0(k,j) - *p * *x;
                    b0(k+1,j) = b0(k+1,j) - *p * *y; 
                }
                if (en <= k+3) {
                    j = en;
                } else {
                    j = k+3;
                }
//c     .......... column modification ..........
                for (i = l; i <= j; i++) {
                    *p = *x * b0(i,k) + *y * b0(i,k+1);
                    b0(i,k) = b0(i,k) - *p;
                    b0(i,k+1) = b0(i,k+1) - *p * *q;
                }
            }
        }

}
