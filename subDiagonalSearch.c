#include<math.h>
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define b0(i,j) B[(i) + (j) * n]
/**
 *  This file is responsible for the sub-diagonal search behavior. 
 *  The intention is that it will be a translation of 
 *  label 70 to right before label 100
 */
/**
 * Only s will be changed after execution.
 * The return value, l, is either low or the
 * index where, for a block like [a b\\c d] we
 * have c being small enough such that a*d + c = a*d
 * in machine arithmetic
 */
int subDiagonalSearch(int n, int low, double* B, int en, double norm, double* s)
{
    int l;
    for (int ll = low; ll <= en; ll++) {
        l = en + low - ll;
        if (l == low) break;
        *s = fabs(b1(l-1,l-1)) + fabs(b1(l,l));
        if (*s == 0.0) *s = norm;
        double tst1 = *s;
        double tst2 = tst1 + fabs(b1(l,l-1));
        if (tst1 == tst2) break;
    }
    return l;
}
