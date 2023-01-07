#include<math.h>
#define a0(i,j) A[(i) + (j) * n]
#define schurMatrix0(i,j) schurMatrix[(i) + (j) * n]
extern int qrIteration(int n, double* A, int en, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m, int schurVectorFlag, int low, int igh, double* schurMatrix);

extern int formShift(int n, int low, double* A, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w);

extern int subDiagonalSearch(int n, int low, double* A, int en, double norm, double* s);

extern int doubleSubDiagonalSearch(int n, double* A, int en, int l, double* s, double x,
        double y, double w, double* p, double* q, double* r, double* zz);

extern void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci); 
