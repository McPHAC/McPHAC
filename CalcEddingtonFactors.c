#include "includes.h"

// Calculate the Eddington factors from P

void CalcEddingtonFactors(column_type  *c, double *mu, double *dmu, int k, int ndepths, int nmu) {
	int i, j;

	for (i=0; i<ndepths; i++) {
		c[i].f[k] = 0.0;
		c[i].h[k] = 0.0;
		for (j=0;j<nmu;j++) {
			c[i].f[k] += mu[j]*mu[j]*c[i].u[j][k]*dmu[j];
			c[i].h[k] += mu[j]*c[i].u[j][k]*dmu[j];
		}
		c[i].f[k] /= c[i].J[k];
		c[i].h[k] /= c[i].J[k];
		for (j=0;j<nmu;j++) {
			c[i].JanisT[j][k] = 3.0/8.0 * ((3.0-mu[j]*mu[j])+(3.0*mu[j]*mu[j]-1)*c[i].f[k]) * c[i].J[k];
			c[i].q[k] += dmu[j]*c[i].JanisT[j][k];
		}
		c[i].q[k] /= c[i].J[k];
	}
}
