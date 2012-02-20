#include "includes.h"

void CalcMeanIntensityJ(column_type  *c, double *dmu, int k, int ndepths, int nmu) {
	// Calculates the mean intensity J at each depth

	int  i, j;

	for (i=0; i<ndepths; i++) {
		c[i].J[k]=0;
		for (j=0; j<nmu; j++) {
			c[i].J[k] += dmu[j]*c[i].u[j][k];
		}
	}
}


