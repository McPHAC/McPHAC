#include "includes.h"

void
CalcCF(double **CF, int i, int k, column_type *c, int ndepths) {
	// Calculates the elements of the C matrix in the Feautrier solution to the radiative transfer

	extern double  *mu;
	int      j, jj;
	double prefactor, delta1, delta2, pi, C, cval, scatterm, cross;

	for (j=0; j<NMU; j++) {
		memset(CF[j], 0.0, NMU*sizeof(double));
	}

	for (j=0; j<NMU; j++) {
		if (i==0) {
			if (AUER) {
				delta2 = c[i].dtau[k]/fabs(mu[j]);
				cval = -1.0/delta2;
				C = delta2/3.0/2.0;
				for (jj=0; jj<NMU; jj++) {
					cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
					scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
					CF[j][jj] =  -c[i+1].rho_opacity[k]*scatterm * C;
				}
				CF[j][j] +=  cval + C;
			} else {
				delta2 = (c[i+1].k[k] + c[i].k[k]) * (c[i+1].y - c[i].y);
				prefactor = 8.0 * mu[j]*mu[j];
				CF[j][j] =  -prefactor/delta2/delta2;
			}
		} else if (i==(ndepths-1)) {
			CF[j][j] = 0.0;
		} else {
			if (AUER) {
				delta1 = c[i-1].dtau[k]/fabs(mu[j]);
				delta2 = c[i].dtau[k]/fabs(mu[j]);
				pi = 0.5*(delta1 + delta2);
				cval = -1.0/pi/delta2;
				C = 1.0/6.0*(1.0+0.5*cval*delta1*delta1);
				for (jj=0; jj<NMU; jj++) {
					cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
					scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
					CF[j][jj] =  -c[i+1].rho_opacity[k]*scatterm * C;
				}
				CF[j][j] +=  cval + C;
			} else {
				delta1 = (c[i].k[k] + c[i-1].k[k]) * (c[i].y - c[i-1].y);
				delta2 = (c[i+1].k[k] + c[i].k[k]) * (c[i+1].y - c[i].y);
				pi = delta1 + delta2;
				prefactor = 8.0 * mu[j]*mu[j];
				CF[j][j] =  -prefactor/pi/delta2;
			}
		}
	}
}
