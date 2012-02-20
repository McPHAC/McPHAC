#include "includes.h"

void
CalcBF(double **BF, int i, int k, column_type *c, int ndepths) {
	// Calculates the elements of the B matrix in the Feautrier solution to the radiative transfer

	extern double  *mu;
	int j, jj;
	double prefactor, delta1, delta2, pi, A, C, a, cval, scatterm, cross;

	for (j=0; j<NMU; j++) {
		memset(BF[j], 0.0, NMU*sizeof(double));
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
					BF[j][jj] =  -c[i].rho_opacity[k]*scatterm * 2.0*C;
				}
				BF[j][j] +=  -cval + 1.0 + 2.0*C;
			} else {
				delta2 = (c[i+1].k[k] + c[i].k[k]) * (c[i+1].y - c[i].y);
				prefactor = 8.0 * mu[j]*mu[j];
				C = -prefactor/delta2/delta2;
				BF[j][j] +=  1.0 - C + 4*fabs(mu[j])/delta2;
			}
		} else if (i==(ndepths-1)) {
			if (AUER) {
				if (DIFFUSION) {
					BF[j][j] =  1.0;
				} else {
					delta1 = c[i-1].dtau[k]/fabs(mu[j]);
					a = -1.0/delta1;
					A = delta1/3.0/2.0;
					for (jj=0; jj<NMU; jj++) {
						cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
						scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
						BF[j][jj] =  -c[i].rho_opacity[k]*scatterm * 2.0*A;
					}
					BF[j][j] +=  -a + 1.0 + 2.0*A;
				}
			} else {
				prefactor = 2.0*fabs(mu[j])/(c[i].y-c[i-1].y);
				BF[j][j] =  1.0 + prefactor/(c[i].k[k] + c[i-1].k[k]);
			}
		} else {
			if (AUER) {
				delta1 = c[i-1].dtau[k]/fabs(mu[j]);
				delta2 = c[i].dtau[k]/fabs(mu[j]);
				pi = 0.5*(delta1 + delta2);
				a = -1.0/pi/delta1;
				cval = -1.0/pi/delta2;
				A = 1.0/6.0*(1.0+0.5*a*delta2*delta2);
				C = 1.0/6.0*(1.0+0.5*cval*delta1*delta1);
				for (jj=0; jj<NMU; jj++) {
					cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
					scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
					BF[j][jj] =  -c[i].rho_opacity[k]*scatterm * (1.0 - A - C);
				}
				BF[j][j] +=  -a - cval + (1.0 - A - C);
			} else {
				delta1 = (c[i].k[k] + c[i-1].k[k]) * (c[i].y - c[i-1].y);
				delta2 = (c[i+1].k[k] + c[i].k[k]) * (c[i+1].y - c[i].y);
				pi = delta1 + delta2;
				prefactor = 8.0 * mu[j]*mu[j];
				A = -prefactor/pi/delta1;
				C = -prefactor/pi/delta2;
				for (jj=0; jj<NMU; jj++) {
					cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
					scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
					BF[j][jj] =  -c[i].rho_opacity[k]*scatterm;
				}
				BF[j][j] +=  1.0 - A - C;
			}
		}
	}
}



