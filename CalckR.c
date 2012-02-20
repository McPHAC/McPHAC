#include "includes.h"

void
CalckR(column_type  *c) {
    // Calculates the Rosseland mean opacity

    double sum, denom, denom2;
    double T;
    int    k;
    void  IntegrateRosselandOpacity();

    sum=0;

    T = pow(10.0,c->logT);
    sum = 0.0;
    denom = 0.0;
    for (k=0; k<NFREQ; k++) {
	sum += dBdT(nu[k],T)/(c->k[k])*dnu[k];
	denom += dBdT(nu[k],T)*dnu[k];
    }
    denom2 = (4.0*SIGMA*pow(T,3.0)/M_PI);
    if (fabs(1.0-denom/denom2) > 0.01) {
	fprintf(stderr, "\nRosseland Mean opacity may be inaccurate at column depth y=%e: denom=%e != %e=4*sigma*T^3/pi\n", c->y, denom, denom2);
    }
    c->kR = denom/sum;
}
