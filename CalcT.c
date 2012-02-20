#include "includes.h"

void CalcT(double **T, int j, int k, column_type *c, int ndepths) {
	// Calculates the elements of the T matrix in the Rybicki solution to the radiative transfer

	extern double  *mu;
	int     icolumn;
	double  prefactor, A, B, C, delta1, delta2, pi, kprefactor;

	for (icolumn=0; icolumn<ndepths; icolumn++) {
		memset(T[icolumn], 0, ndepths*sizeof(double));

		if (icolumn==0) { // Top of the atmosphere
			prefactor= fabs(mu[j])/(c[icolumn+1].y-c[icolumn].y);
			A = 0.0;
			C =  -prefactor/(0.5*(c[icolumn].k[k] + c[icolumn+1].k[k]));
			B =  1.0 - C;
		}  else if (icolumn<ndepths-1) { // Inside the atmosphere
			prefactor=8.0 * mu[j]*mu[j];
			delta1 = (c[icolumn].k[k] + c[icolumn-1].k[k]) * (c[icolumn].y - c[icolumn-1].y);
			delta2 = (c[icolumn+1].k[k] + c[icolumn].k[k]) * (c[icolumn+1].y - c[icolumn].y);
			pi = delta1 + delta2;
			A =  - prefactor/pi/delta1;
			C =  - prefactor/pi/delta2;
			B = 1.0 -  A -  C;
		} else { // Bottom of the atmosphere
			prefactor= 2.0*fabs(mu[j])/(c[icolumn].y-c[icolumn-1].y);
			kprefactor = prefactor/(c[icolumn].k[k] + c[icolumn-1].k[k]);
			C=0;
			A= - kprefactor;
			B = 1.0 - A;
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


