#include "includes.h"

void CalcTeff(column_type  *c, double *nu, double *dnu, int ndepths, int nfreq) {
    // Calculates the effective temperature at each depth from the flux

    int i, k;

    for (i=0; i<ndepths; i++) {
        c[i].Teff = 0.0;
	for (k=0; k<nfreq; k++) {
            c[i].Teff += c[i].F[k]*dnu[k];
	}
        c[i].Teff = pow(c[i].Teff/SIGMA,0.25);
    }
}
