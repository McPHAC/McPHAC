#include "includes.h"

void
UpdateW(double **W, double** V, double** TInverse, double** U, int ndepths) {
    // Calculates the W matrix in the Rybicki method

    int  i, s;
    double   **tmp1, **tmp2;

    tmp1 = ddvector(ndepths);
    tmp2 = ddvector(ndepths);
    for (i=0; i<ndepths; i++) {
        tmp1[i]=dvector(ndepths);
        tmp2[i]=dvector(ndepths);
    }

    MatrixDiagonalMatrixMultiply(tmp2, TInverse, U, ndepths);
    DiagonalMatrixMatrixMultiply(tmp1,V,tmp2, ndepths);

    for (i=0; i<ndepths; i++) {
        for (s=0; s<ndepths; s++) {
            W[i][s] -= tmp1[i][s];
        }
    }

    for (i=0; i<ndepths; i++) {
        free(tmp1[i]);
        free(tmp2[i]);
    }
    free(tmp1);
    free(tmp2);
}

