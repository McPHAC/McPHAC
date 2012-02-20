#include "includes.h"

void DiagonalMatrixMatrixMultiply(double** ans, double** A, double** B, int n) {
    // Multiplies diagonal matrix A with matrix B, both of dimension n by n
    int  i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    ans[i][j] = A[i][i] * B[i][j];
	}
    }
}
