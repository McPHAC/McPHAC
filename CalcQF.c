#include "includes.h"

void
CalcQF(double *QF, int i, int k, int knu, column_type *c, int ndepths) {
    // Calculates the elements of the Q vector in the Feautrier solution to the radiative transfer
    extern double  *nu, *mu;
    int      j;
    double   T, T1, T2, delta2, Q, Q1, Q2, A, B, C, a, cval, delta1, pi;

    memset(QF, 0.0, NMU*sizeof(double));

    T = pow(10.0, c[i].logT);
    for (j=0; j<NMU; j++) {
        if (i==0) {
            T2 = pow(10.0, c[i+1].logT);
            if (AUER) {
		delta2 = c[i].dtau[knu]/fabs(mu[j]);
                C = delta2/3.0/2.0;
                B = 2.0*C;
                Q = (1.0 - c[i].rho_opacity[knu]) * Bnu(nu[k],T);
                Q2 = (1.0 - c[i+1].rho_opacity[knu]) * Bnu(nu[k],T2);
                QF[j] = B*Q1 + C*Q2;
            } else {
                QF[j] = 0.0;
            }
        } else if (i==(ndepths-1)) {
	    T1 = pow(10.0, c[i-1].logT);
            if (AUER) {
                if (DIFFUSION) {
		    QF[j] = Bnu(nu[k],T);
                } else {
		    delta1 = c[i-1].dtau[knu]/fabs(mu[j]);
		    A = delta1/3.0/2.0;
		    B = 2.0 * A;
		    Q = (1.0 - c[i].rho_opacity[knu]) * Bnu(nu[k],T);
		    Q1 = (1.0 - c[i-1].rho_opacity[knu]) * Bnu(nu[k],T1);
		    QF[j] =  A*Q1 + B*Q;
		    QF[j] += Bnu(nu[k],T);
                }
            } else {
		QF[j] = Bnu(nu[k],T);
            }
        } else {
	    T1 = pow(10.0, c[i-1].logT);
	    T2 = pow(10.0, c[i+1].logT);
            if (AUER) {
		delta1 = c[i-1].dtau[knu]/fabs(mu[j]);
		delta2 = c[i].dtau[knu]/fabs(mu[j]);
                pi = 0.5*(delta1 + delta2);
                a = 1.0/pi/delta1;
                cval = 1.0/pi/delta2;
                A = 1.0/6.0*(1.0+0.5*a*delta2*delta2);
		A *= 0.0;
                C = 1.0/6.0*(1.0+0.5*cval*delta1*delta1);
		C *= 0.0;
                B = 1.0 - A - C;
                Q = (1.0 - c[i].rho_opacity[knu]) * Bnu(nu[k],T);
                Q1 = (1.0 - c[i-1].rho_opacity[knu]) * Bnu(nu[k],T1);
                Q2 = (1.0 - c[i+1].rho_opacity[knu]) * Bnu(nu[k],T2);
                QF[j] =  A*Q1 + B*Q + C*Q2;
            } else {
                QF[j] = (1.0 - c[i].rho_opacity[knu]) * Bnu(nu[k],T);
            }
        }
    }
}
