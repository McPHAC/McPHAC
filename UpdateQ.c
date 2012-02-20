#include "includes.h"

void UpdateQ(double *Q, double **V, double **TInverse, double *K, int ndepths) {
    // Calculates the vector Q in the Rybicki method

    int  i;
    double   *tmp2, *tmp3;

    tmp2 = dvector(ndepths);
    tmp3 = dvector(ndepths);

    MatrixVectorMultiply(tmp2, TInverse, K, ndepths);
    MatrixVectorMultiply(tmp3,V,tmp2, ndepths);

    for (i=0; i<ndepths; i++) {
        Q[i] -= tmp3[i];
    }

    free(tmp3);
    free(tmp2);
}
