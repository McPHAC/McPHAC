#include "includes.h"

void
CalcKt(double *K, int j, int k, column_type *c, int ndepths) {
    // Calculates the elements of the K vector in the Rybicki solution in the temperature correction proccedure

    extern double  *nu;
    int      i, l;
    double   T;
    double   num, denom;


    for (i=0; i<(ndepths); i++) {
        if (i==0) {
            K[i] = 0.0;
        } else if (i==(ndepths-1)) {
            T = pow(10.0, c[i].logT);
            K[i] = Bnu(nu[k],T);
        } else {
            T = pow(10.0, c[i].logT);
            num = 0.0;
            denom = 0.0;
            for (l=0; l<NFREQ; l++) {
                num += Bnu(nu[l], T)*c[i].kappa[l]*dnu[l];
                denom += dBdT(nu[l], T)*c[i].kappa[l]*dnu[l];
            }
            K[i] = (Bnu(nu[k],T) - dBdT(nu[k],T)*num/denom)*(1.0-c[i].rho_opacity[k]);
        }
    }
}



