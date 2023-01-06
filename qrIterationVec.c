#define b0(i,j) B[(i) + (j) * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define z0(i,j) Z[(i) + (j) * n]
#define z1(i,j) Z[(i - 1) + (j - 1) * n]
#include<stdbool.h>
#include<math.h>
/**
 * n: One dimension 
 */
void qrIterationVec(int n, double* B, int en, int l, double* s, double* x, 
        double * y, double* p, double* q, double* r, double* zz, int m, int low,
        int igh, double *Z)
{
    int i,j,k;
    bool notLast;
    for (k = m; k <= en - 1; k++) {
        notLast = k != en - 1;
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
//c     .......... row modification ..........
            for (j = k; j < n; j++) {
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
            for (i = 0; i <= j; i++) {
                *p = *x * b0(i,k) + *y * b0(i,k+1) + *zz * b0(i,k+2);
                b0(i,k) = b0(i,k) - *p;
                b0(i,k+1) = b0(i,k+1) - *p * *q;
                b0(i,k+2) = b0(i,k+2) - *p * *r;
            }
//c     .......... accumulate transformations  ..........
            for (i = low; i <= igh; i++) {
                *p = *x * z0(i,k) + *y * z0(i,k+1) + *zz * z0(i,k+2);
                z0(i,k) = z0(i,k) - *p;
                z0(i, k + 1) = z0(i,k + 1) - *p * *q;
                z0(i, k + 2) = z0(i,k + 2) - *p * *r;
            }
        } else {
//c     .......... row modification ..........
            for (j = k; j < n; j++) {
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
            for (i = 0; i <= j; i++) {
                *p = *x * b0(i,k) + *y * b0(i,k+1);
                b0(i,k) = b0(i,k) - *p;
                b0(i,k+1) = b0(i,k+1) - *p * *q;
            }
//c     .......... accumulate transformations  ..........
            for (i = low; i <= igh; i++) {
                *p = *x * z0(i,k) + *y * z0(i,k+1);
                z0(i,k) = z0(i,k) - *p;
                z0(i, k + 1) = z0(i,k+1) - *p * *q;
            }
        }
    }
}
