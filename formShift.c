#define b0(i,j) B[(i) + (j) * n]
#define b1(i,j) B[(i - 1) + (j - 1) * n]
#include<math.h>
/**
 *  This file is responsible for the form shift behavior. 
 *  The intention is that it will be a translation of 
 *  label 100 to right before label 130 
 */
/**
 *  Inputs that are passed as a pointer are modified after calling
 *  this function. 
 *  Inputs:
 *      n - Size of the matrix B
 *      low - The lower bound of 
 *  Outputs:
 */
int formShift(int n, int low, double* B, int* ierr, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w)
{
    int retVal = 0;
    int ione = 1;
    int na = en - 1;
    int enm2 = na - 1;
    *x = b1(en,en);
    // We have found a single root, so we return 1.
    if (l == en) return 1;
    *y = b1(na,na);
    *w = b1(en,na) * b1(na,en);
    // This means we have found a double root, so we return 2
    if (l == na) return 2;
    // This means we errored out, so we set the error index variable, then return a 3
    if (itn == 0) {
        *ierr = en;
        return 3;
    }
    // if its is 10 or 20, then we perform an "exceptional" shift, otherwise
    // we terminate successfully
    if (its != 10 && its != 20) return 0;
    *t = *t + *x;
    for (int i = low; i <= en; i++) 
        b1(i,i) -= *x;
    *s = fabs(b1(en,na)) + fabs(b1(na,enm2));
    *x = 0.75 * *s;
    *y = *x;
    *w = -0.4375 * *s * *s;
    return 0;
}
