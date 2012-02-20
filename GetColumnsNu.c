#include "includes.h"

column_type* GetColumnsNu(column_type *c, int ndepths, int nmunu, int k, int knu, int *ndepthsnu)
{
    // Calculates the frequency specific column_type array based on the colun_type array c
    int i, s;
    column_type *cnu;
    double  *y, *ynu, *T, *Tnu, dy, dlogy;

    cnu = AllocateColumns(*ndepthsnu, nmunu, 1);
    y = (double *) malloc((ndepths)*sizeof(double));
    T = (double *) malloc((ndepths)*sizeof(double));
    for (i=0; i<ndepths; i++) { // Copy values from c into arrays
        y[i] = pow(10.0, c[i].logy);
        T[i] = pow(10.0, c[i].logT);
    }
    ynu = (double *) malloc((*ndepthsnu)*sizeof(double));
    Tnu = (double *) malloc((*ndepthsnu)*sizeof(double));
    for (s=0; (s<ndepths-1)?(c[s].tau[k]<MAXTAU):0; s++); // Find first depth with tau>MAXTAU
    dy = (y[s]-y[0])/(double)(*ndepthsnu-1);
    dlogy = 1.0*(c[s].logy-c[0].logy)/(1.0*(*ndepthsnu-1));
    for (i=0; i<*ndepthsnu; i++) { // Set column depths to use in cnu
        if (USELOGCOLNU) {
	    ynu[i] = pow(10.0, c[0].logy + dlogy*i); // Logarithmically spaced
        } else {
            ynu[i] = y[0] + dy*i; // Linearly spaced
        }
    }
    InterpolateArray(y, T, ndepths, ynu, Tnu, *ndepthsnu, 0, 0.0, 0, 0.0);

    for (i=0; i<*ndepthsnu; i++) { // Set column depths to use in cnu
        cnu[i].y = ynu[i];
        cnu[i].logy = log10(cnu[i].y);
        cnu[i].logT = log10(Tnu[i]);
        cnu[i].pressure = GSURFACE * cnu[i].y;
        cnu[i].rho = CalcRho(cnu[i].pressure, pow(10.0,cnu[i].logT));
        CalcOpacities(&cnu[i], knu, 1, k);
    }
    cnu[0].tau[knu] = cnu[0].k[knu]*cnu[0].y;
    for (i=0; i+1<*ndepthsnu; i++) {
        cnu[i].dy = cnu[i+1].y-cnu[i].y;
        cnu[i].dtau[knu] = 0.5*(cnu[i].k[knu]+cnu[i+1].k[knu])*cnu[i].dy;
        cnu[i+1].tau[knu] = cnu[i].tau[knu] + cnu[i].dtau[knu];
    }
    cnu[*ndepthsnu-1].dy = 0.0;
    cnu[*ndepthsnu-1].dtau[knu] = 0.0;

    free(y);
    free(T);
    free(ynu);
    free(Tnu);

    return (cnu);
}
