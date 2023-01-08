#define b0(i,j) B[(i) + (j) * n]
#include<math.h>
/**
 *  Inputs that are passed as a pointer are modified after calling
 *  this function. 
 *  the index m is returned after this 
 *
 *  Inputs:
 *  n: The dimension of our matrix B
 *  B: The Matrix we are searching for small sub-diagonal elements
 *  en: The index of the eigenvalue we are currently trying to find
 *  l: The index that was returned by subDiagonalSearch
 *  s: Value that is used for testing for small elements in other functions
 *  x: Value used to help construct our testing criteria
 *  y: Value used to help construct our testing criteria
 *  w: Value used to help construct our testing criteria
 *  p: Value used to help construct our testing criteria
 *  q: Value used to help construct our testing criteria
 *  r: Value used to help construct our testing criteria
 *  zz: Value used to help construct our testing criteria
 *
 *  Return:
 *  m: The index under which we have as a candidate for a double root
 */

int doubleSubDiagonalSearch(int n, double* B, int en, int l, double* s, 
        double x, double y, double w, double* p, double* q, double* r, double* zz
        )
{
    int m;
    for (m = en - 2; m >= l; m--) {
        *zz = b0(m,m);
        *r = x - *zz;
        *s = y - *zz;
        *p = (*r * *s - w) / b0(m+1,m) + b0(m,m+1);
        *q = b0(m+1,m+1) - *zz - *r - *s;
        *r = b0(m+2,m+1);
        *s = fabs(*p) + fabs(*q) + fabs(*r);
        *p = *p / *s;
        *q = *q / *s;
        *r = *r / *s;
        if (m == l) break;
        double tst1 = fabs(*p) * (fabs(b0(m-1, m-1)) + fabs(*zz) + fabs(b0(m+1,m+1)));
        double tst2 = tst1 + fabs(b0(m,m-1)) * (fabs(*q) + fabs(*r));
        if (tst1 == tst2) break;
    }
    for (int i = m + 2; i <= en; i++) {
        b0(i,i-2) = 0.0;
        if (i == m + 2) continue;
        b0(i,i-3) = 0.0;
    }
    return m;
}
