#include "includes.h"

void CalcJtbar(double **W, double *Q, column_type  *c, int ndepths) {
    // Calculates the quantity Jbar, used in the Rybicki solution of the temperature correction proceedure
    extern double  *dnu;

    int  i, k;
    double **tmp1, *tmp2, nutot;

    tmp1 = ddvector(ndepths);
    tmp2 = dvector(ndepths);
    for (i=0; i<ndepths; i++) {
        tmp1[i]=dvector(ndepths);
    }

    InvertMatrix(tmp1, W, ndepths);
    MatrixVectorMultiply(tmp2, tmp1, Q, ndepths);

    nutot=0;
    for (k=0; k<NFREQ; k++) {
        nutot+=dnu[k];
    }

    for (i=0; i<ndepths; i++) {
        c[i].Jtbar = tmp2[i];
    }

    for (i=0; i<ndepths; i++) {
        free(tmp1[i]);
    }
    free(tmp1);
    free(tmp2);

 }
