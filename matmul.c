#define a0(i,j) A[(i) + (j) * nA]
#define b0(i,j) B[(i) + (j) * nB]
#define c0(i,j) C[(i) + (j) * nB]
#include<stdlib.h>
/**
 * This function takes in two matrices and returns the 
 * product A * B. We need to know the dimensions thereof. If
 * the user gives incorrect dimensions, then the result will
 * contain garbage data.
 *
 * This is a general purpose algorithm for dense matrices.
 *
 * Inputs:
 * A - nA by mA matrix of double precision elements.
 * nA - number of rows for the matrix A
 * mA - number of cols for the matrix A
 * B - nB by mB matrix of double precision elements.
 * nB - number of rows for the matrix B
 * mB - number of cols for the matrix B
 *
 * Outputs:
 * C - nA by mB matrix of double precision elements such that
 * C = A * B, where * is the standard matrix multiplication.
 * if mA is not equal to nB, then matrix multiplication
 * is not defined, so we return null
 * 
 * Note: we use malloc to create C, so if it is not null, 
 * the user must remember to free the memory.
 */
double *matmul(double *A, int nA, int mA, double *B, int nB, int mB) {
    if (mA != nB) {
        // This means multiplication is not defined, so return null
        return NULL;
    }
    double *C = (double *) malloc(nA * mB * sizeof(double)); 
    for (int i = 0; i < nA; i++) {
        for (int j = 0; j < mB; j++) {
            double sum = 0.0;
            for (int k = 0; k < mA; k++) {
                sum += a0(i,k) * b0(k,j);
            }
            c0(i,j) = sum;
        }
    }
    return C;
}
