#include "includes.h"

void
CalcK(double *K, int j, int k, int knu, column_type *c, int ndepths) {
    // Calculates the elements of the K vector in the Rybicki solution to the radiative transfer

    extern double  *nu;
    int      i;
    double   T;
    double   dT;
    double   dBdtau;


    for (i=0; i<(ndepths); i++) {
        if (i==0) {
            K[i] = 0.0;
        } else if (i==(ndepths-1)) {
            T = pow(10.0, c[i].logT);
            dT = pow(10.0, c[i-1].logT) - T;
            dBdtau = (Bnu(nu[k],T) - Bnu(nu[k],T+dT)) / 0.5 / (c[i].k[knu] + c[i-1].k[knu]) / (c[i].y - c[i-1].y);
            K[i] = Bnu(nu[k],T);
        } else {
            T = pow(10.0, c[i].logT);
            K[i] = (1.0 - c[i].rho_opacity[knu]) * Bnu(nu[k],T);
        }
    }
}
