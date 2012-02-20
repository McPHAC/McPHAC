#include "includes.h"

void CalcTt(double **T, int j, int k, column_type *c, int ndepths) {
    // Calculates the elements of the T matrix in the Rybicki solution in the temperature correction proccedure
    int     icolumn;
    double  prefactor, A, B, C, delta1, delta2, pi;

    for (icolumn=0; icolumn<ndepths; icolumn++) {
        memset(T[icolumn], 0, ndepths*sizeof(double));

        if (icolumn==0) { // Top of the atmopshere
            prefactor= 1.0/(c[icolumn+1].y-c[icolumn].y);
            A = 0.0;
            C =  -prefactor*c[icolumn+1].f[k]/(0.5*(c[icolumn].k[k] + c[icolumn+1].k[k]));
            B =  c[icolumn].h[k] + prefactor*c[icolumn].f[k]/(0.5*(c[icolumn].k[k] + c[icolumn+1].k[k]));
        }  else if (icolumn<ndepths-1) { // Inside the atmosphere
            prefactor = 8.0;
            delta1 = (c[icolumn].k[k] + c[icolumn-1].k[k]) * (c[icolumn].y - c[icolumn-1].y);
            delta2 = (c[icolumn+1].k[k] + c[icolumn].k[k]) * (c[icolumn+1].y - c[icolumn].y);
            pi = delta1 + delta2;
            A =  - c[icolumn-1].f[k]*prefactor/pi/delta1;
            C =  - c[icolumn+1].f[k]*prefactor/pi/delta2;
            B = 1.0*(1.0-c[icolumn].rho_opacity[k]) + c[icolumn].f[k]*prefactor/pi/delta1 + c[icolumn].f[k]*prefactor/pi/delta2;
        } else { // Bottom of the atmosphere
            C = 0.0;
            A = 0.0;
            B = 1.0;
        }


        if (icolumn>0) {
            T[icolumn][icolumn-1]=A;
        }

        T[icolumn][icolumn]=B;

        if (icolumn<(ndepths-1)) {
            T[icolumn][icolumn+1]=C;
        }
    }
}



