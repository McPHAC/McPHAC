#include "includes.h"

void CalcFluxes(column_type  *c, double *mu, double *dmu, int k, int ndepths, int nmu) {
    // Calculates the flux through each depth

    int i,j;

    c[0].F[k] = 0.0;
    for (j=0; j<nmu; j++) {
        c[0].F[k] += M_PI*4.0*c[0].u[j][k]*fabs(mu[j])*dmu[j];
    }
    for (i=1; i<ndepths-1; i++) {
        c[i].F[k] = 4.0*M_PI*(c[i+1].f[k]*c[i+1].J[k]-c[i].f[k]*c[i].J[k])/c[i].dtau[k];
    }
    c[ndepths-1].F[k] = 0.0;
}
