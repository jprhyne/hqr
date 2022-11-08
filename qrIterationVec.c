#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define z1(i,j) Z[(i - 1) + (j - 1) * n]
#include<stdbool.h>
#include<math.h>
extern void hqr_qriter_(int* nm, int* n, double* h, int* en, int* na,
        int* l, double* s, double* x, double* y,
        double* p, double* q, double* r, double* zz, int* m);
/**
 * n: One dimension 
 */
void qrIterationVec(int n, double* B, int en, int na, int l, double* s, double* x, 
        double * y, double* p, double* q, double* r, double* zz, int m, int low,
        int igh, double *Z)
{
    int i,j,k;
    bool notLast;
    for (k = m; k <= na; k++) {
        notLast = k != na;
        if (k != m) {
            *p = b1(k,k-1);
            *q = b1(k+1,k-1);
            *r = 0.0;
            if (notLast) {
                *r = b1(k+2,k-1);
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
                b1(k,k-1) = -*s * *x;
            } else if (l != m) {
                b1(k,k-1) = -b1(k,k-1);
            }
            *p = *p + *s;
            *x = *p / *s;
            *y = *q / *s;
            *zz = *r / *s;
            *q = *q / *p;
            *r = *r / *p;
            if (notLast) {
//c     .......... row modification ..........
                for (j = k; j <= en; j++) {
                    *p = b1(k,j) + *q * b1(k+1,j) + *r * b1(k+2,j);
                    b1(k,j) = b1(k,j) - *p * *x;
                    b1(k+1,j) = b1(k+1,j) - *p * *y;
                    b1(k+2,j) = b1(k+2,j) - *p * *zz;
                }
                if (en <= k+3) {
                    j = en;
                } else {
                    j = k+3;
                }
//c     .......... column modification ..........
                for (i = l; i <= j; i++) {
                    *p = *x * b1(i,k) + *y * b1(i,k+1) + *zz * b1(i,k+2);
                    b1(i,k) = b1(i,k) - *p;
                    b1(i,k+1) = b1(i,k+1) - *p * *q;
                    b1(i,k+2) = b1(i,k+2) - *p * *r;
                }
//c     .......... accumulate transformations  ..........
                for (i = low; i <= igh; i++) {
                    *p = *x * z1(i,k) + *y * z1(i,k+1);
                    z1(i,k) = z1(i,k) - *p;
                    z1(i, k + 1) = z1(1,k+1) - *p * *q;
                }
            } else {
//c     .......... row modification ..........
                for (j = k; j <= en; j++) {
                    *p = b1(k,j) + *q * b1(k+1,j);
                    b1(k,j) = b1(k,j) - *p * *x;
                    b1(k+1,j) = b1(k+1,j) - *p * *y; 
                }
                if (en <= k+3) {
                    j = en;
                } else {
                    j = k+3;
                }
//c     .......... column modification ..........
                for (i = l; i <= j; i++) {
                    *p = *x * b1(i,k) + *y * b1(i,k+1);
                    b1(i,k) = b1(i,k) - *p;
                    b1(i,k+1) = b1(i,k+1) - *p * *q;
                }
//c     .......... accumulate transformations  ..........
                for (i = low; i <= igh; i++) {
                    *p = *x * z1(i,k) + *y * z1(i,k+1) + *zz * z1(i,k+2);
                    z1(i,k) = z1(i,k) - *p;
                    z1(i, k + 1) = z1(1,k+1) - *p * *q;
                    z1(i, k + 1) = z1(1,k+1) - *p * *r;
                }
            }
        }

}