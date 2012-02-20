#include "includes.h"

void CalcV(double **V, int j, int k, column_type *c, int ndepths) {
    // Calculates the elements of the V matrix in the Rybicki solution to the radiative transfer

    extern double  *dmu;
    double   weight;
    int        i;

    weight = dmu[j];
    for (i=0; i<ndepths; i++) {
        memset(V[i], 0, ndepths*sizeof(double));
        V[i][i] = weight;
    }
}

