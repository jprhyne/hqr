#include<math.h>
#define a1(i,j) A[(i - 1) + (j - 1) * n]
#define schurMatrix1(i,j) schurMatrix[(i - 1) + (j - 1) * n]
extern int qrIteration(int n, double* A, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m);

extern int qrIterationVec(int n, double* A, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m, int low, int igh, double* schurMatrix);

extern int formShift(int n, int low, double* A, int* ierr, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w);

extern int subDiagonalSearch(int n, int low, double* A, int en, double norm, double* s);

extern int doubleSubDiagonalSearch(int n, double* A, int en, int enm2, int l, double* s, double x,
        double y, double w, double* p, double* q, double* r, double* zz);

extern void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci); 
