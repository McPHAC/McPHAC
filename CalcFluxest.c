#include "includes.h"

void CalcFluxest(column_type  *c, double *mu, double *dmu, int ndepths, int nfreq) {
    // Calculates flux through each depth based on the temperature correction J (Jt)

    int i,k;

    for (k=0; k<nfreq; k++) {
	for (i=0; (i+1)<ndepths; i++) {
	    c[i].F[k] = 4.0*M_PI*(c[i+1].f[k]*c[i+1].Jt[k]-c[i].f[k]*c[i].Jt[k])/c[i].dtau[k];
	}
	c[ndepths-1].F[k] = c[ndepths-2].F[k];
    }
}
