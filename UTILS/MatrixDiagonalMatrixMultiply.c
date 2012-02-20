#include "includes.h"

void MatrixDiagonalMatrixMultiply(double** ans, double** A, double** B, int n) {
    // Multiplies matrix A with diagonal matrix B, both of dimension n by n
    int  i, j;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    ans[i][j] = A[i][j] * B[j][j];
	}
    }
}
