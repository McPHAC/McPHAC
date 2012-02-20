#include "includes.h"

void CalcJbar(double **W, double *Q, column_type  *c, int ndepths) {
    // Calculates the quantity Jbar, used in the Rybicki solution of the radiative transfer

    extern double  *dnu;

    int  i, k;
    double **tmp1, *tmp2, nutot,a;

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
        c[i].Jbar = tmp2[i];

        c[i].Jbarb=0;
        for (k=0; k<NFREQ; k++) {
            c[i].Jbarb += c[i].J[k]*dnu[k];
        }
        c[i].Jbarb /= nutot;
        a=c[i].Jbarb/c[i].Jbar;
    }

    for (i=0; i<ndepths; i++) {
        free(tmp1[i]);
    }
    free(tmp1);
    free(tmp2);

}


