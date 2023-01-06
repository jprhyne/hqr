#include<math.h>
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#define b0(i,j) B[(i) + (j) * n]
/**
 *  This file is responsible for the sub-diagonal search behavior. 
 *  The intention is that it will be a translation of 
 *  label 70 to right before label 100
 *
 *  Inputs:
 *  n: The dimension of B
 *  low: The lowest index we are considering for eigenvalues
 *  B: The matrix we are trying to find the eigenvalues of
 *  en: The index of the eigenvalue we are currently trying to find
 *  norm: The sum of the absolute elements of the original matrix A
 *  s: Value used for testing sub-diagonal elements
 *
 *  Return:
 *  l: The index that contains a small sub-diagonal element and thus we
 *      have as a candidate for an eigenvalue
 */
int subDiagonalSearch(int n, int low, double* B, int en, double norm, double* s)
{
    int l;
    for (int ll = low; ll <= en; ll++) {
        l = en + low - ll;
        if (l == low) break;
        *s = fabs(b0(l-1,l-1)) + fabs(b0(l,l));
        if (*s == 0.0) *s = norm;
        double tst1 = *s;
        double tst2 = tst1 + fabs(b0(l,l-1));
        if (tst1 == tst2) break;
    }
    return l;
}
