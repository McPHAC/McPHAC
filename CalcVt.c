#include "includes.h"

void CalcVt(double **V, int j, int k, column_type *c, int ndepths) {
	// Calculates the elements of the V matrix in the Rybicki solution in the temperature correction procedure

	double   weight, denom, T;
	int        i, l;

	for (i=0; i<ndepths; i++) {
		denom = 0.0;
		T = pow(10.0,c[i].logT);
		for (l=0; l<NFREQ; l++) {
			denom += dBdT(nu[l], T)*c[i].kappa[l]*dnu[l];
		}
		weight = c[i].kappa[k]*dnu[k]/denom;
		memset(V[i], 0, ndepths*sizeof(double));
		V[i][i] = weight;
	}
}



