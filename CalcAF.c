#include "includes.h"

void
CalcAF(double **AF, int i, int k, column_type *c, int ndepths) {
    // Calculates the elements of the A matrix in the Feautrier solution to the radiative transfer

    extern double  *mu;
    int      j, jj;
    double prefactor, delta1, delta2, pi, A, a, scatterm, cross;

    for (j=0; j<NMU; j++) {
        memset(AF[j], 0.0, NMU*sizeof(double));
    }

    for (j=0; j<NMU; j++) {
        if (i==0) {
	    AF[j][j] = 0.0;
        } else if (i==(ndepths-1)) {
            if (AUER) {
                if (DIFFUSION) {
          	    AF[j][j] = 0.0;
                } else {
		    delta1 = c[i-1].dtau[k]/fabs(mu[j]);
		    a = -1.0/delta1;
		    A = delta1/3.0/2.0;
		    for (jj=0; jj<NMU; jj++) {
			cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
			scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
			AF[j][jj] =  -c[i-1].rho_opacity[k]*scatterm * A;
		    }
		    AF[j][j] +=  a + A;
                }
	    } else {
		prefactor= 2.0*fabs(mu[j])/(c[i].y-c[i-1].y);
		AF[j][j] = -prefactor/(c[i].k[k] + c[i-1].k[k]);
	    }
        } else {
            if (AUER) {
		delta1 = c[i-1].dtau[k]/fabs(mu[j]);
		delta2 = c[i].dtau[k]/fabs(mu[j]);
		pi = 0.5*(delta1 + delta2);
		a = -1.0/pi/delta1;
		A = 1.0/6.0*(1.0+0.5*a*delta2*delta2);
		for (jj=0; jj<NMU; jj++) {
		    cross = 3.0/8.0*((3.0-mu[j]*mu[j]) + (3.0*mu[j]*mu[j]-1)*mu[jj]*mu[jj]);
		    scatterm = (ANIST) ? cross*dmu[jj] : dmu[jj];
		    AF[j][jj] =  -c[i-1].rho_opacity[k]*scatterm * A;
		}
		AF[j][j] +=  a + A;
            } else {
		delta1 = (c[i].k[k] + c[i-1].k[k]) * (c[i].y - c[i-1].y);
		delta2 = (c[i+1].k[k] + c[i].k[k]) * (c[i+1].y - c[i].y);
		pi = delta1 + delta2;
		prefactor = 8.0 * mu[j]*mu[j];
		AF[j][j] =  -prefactor/pi/delta1;
	    }
        }
    }
}
