#include<math.h>
/**
 *  This function divides two complex numbers in real arithmetic
 *  and stores the result in cr and ci
 *  That is (cr,ci) = (ar,ai)/(br,bi)
 */
void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci) {
    double s,ars,ais,brs,bis;
    s = fabs(br) + fabs(bi);
    ars = ar/s;
    ais = ai/s;
    brs = br/s;
    bis = bi/s;
    s = brs * brs + bis * bis; 
    *cr = (ars*brs + ais*bis) / s;
    *ci = (ais*brs - ars*bis) / s;
}
