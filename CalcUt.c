#include "includes.h"

void CalcUt(double **U,  int k, column_type *c, int ndepths) {
	// Calculates the elements of the U matrix in the Rybicki solution in the temperature correction procedure

	double T;
	int  i;

	memset(U[0], 0, ndepths*sizeof(double));
	U[0][0] = 0.0;   // Top boundary condition
	for (i=1; i<(ndepths-1); i++) {
		memset(U[i], 0, ndepths*sizeof(double));
		T = pow(10.0,c[i].logT);
		U[i][i] = -dBdT(nu[k],T)*(1.0-c[i].rho_opacity[k]);
	}
	memset(U[ndepths-1], 0, ndepths*sizeof(double));
	U[ndepths-1][ndepths-1] = 0.0;   // Bottom boundary condition
}



