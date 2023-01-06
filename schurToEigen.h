#include<math.h>
#define t0(i,j) T[(i) + (j) * n]
#define t1(i,j) T[(i - 1) + (j - 1) * n]
#define eigenMatrix0(i,j) eigenMatrix[(i) + (j) * n]
#define eigenMatrix1(i,j) eigenMatrix[(i - 1) + (j - 1) * n]
extern void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci); 
