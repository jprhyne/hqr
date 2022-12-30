extern int qrIteration(int n, double* B, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m);

extern int qrIterationVec(int n, double* B, int en, int na, int l, double* s,
        double* x, double* y, double* p, double* q, double* r, double* zz,
        int m, int low, int igh, double* eigenMatrix);

extern int formShift(int n, int low, double* B, int* ierr, int its, int itn,
        int en, int l, double* s, double* t, double* x, double* y, double* w);

extern int subDiagonalSearch(int n, int low, double* B, int en, double norm, double* s);

extern int doubleSubDiagonalSearch(int n, double* B, int en, int enm2, int l, double* s, double x,
        double y, double w, double* p, double* q, double* r, double* zz);

extern void hqr2schur_(int *nm, int *n, int *low, int *igh, double *h, double *wr,
        double *wi, double *z, int *ierr);

extern void hqr2eigen_(int *nm, int *n, int *low, int *igh, double *h, double *wr,
        double *wi, double *z, int *ierr);

extern void hqr_(int *nm, int *n, int *low, int *igh, double *h, double *wr,
        double *wi, int *ierr);

extern int hqr(int nm, int n, int low, int igh, double *A, double *eigenValsReal, double *eigenValsImag, int schurVectorFlag, double *eigMat);

extern void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci); 

extern double*matmul(double*A, int nA, int mA, double *B, int nB, int mB);

extern double*matsub(double*A, int nA, int mA, double *B, int nB, int mB);

