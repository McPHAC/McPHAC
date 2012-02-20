#include "includes.h"

void MatrixVectorMultiply(double *ans, double **A, double *V, int n) {
    // Multiplies n by n matrix A by vector V

    int  i, k;
    double  sum;
    for (i=0; i<n; i++) {
        sum=0;
        for (k=0; k<n; k++) {
            sum += A[i][k] * V[k];
        }
        ans[i]=sum;
    }
}
