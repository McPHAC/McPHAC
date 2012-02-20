#include "includes.h"

void InterPolEdd(column_type * cnu, int knu, int ndepthsnu, column_type * c, int k, int ndepths)
// Fills the full column_type array c with Eddington factors based on interpolation of those in cnu
{
    int i,s,n;
    double  *ya, *fa, *ha, *yinterp, *finterp, *hinterp;

    for (s=0; (s<ndepths)?(c[s].y<=cnu[ndepthsnu-1].y):0; s++); // s is index of first depth of c beyond the range of cnu
    n = ndepthsnu + ndepths-s;
    ya = (double *)malloc((n)*sizeof(double));
    fa = (double *)malloc((n)*sizeof(double));
    ha = (double *)malloc((n)*sizeof(double));
    for (i=0; i<ndepthsnu; i++) { // Set y, f, and h values in cnu
	ya[i] = cnu[i].y;
	fa[i] = cnu[i].f[knu];
	ha[i] = cnu[i].h[knu];
    }
    for (i=s; i<ndepths; i++) { // Set y, f, and h values beyond depth s
	ya[ndepthsnu+i-s] = c[i].y;
	fa[ndepthsnu+i-s] = cnu[ndepthsnu-1].f[knu];
	ha[ndepthsnu+i-s] = cnu[ndepthsnu-1].h[knu];
    }
    yinterp = (double *)malloc((ndepths)*sizeof(double));
    finterp = (double *)malloc((ndepths)*sizeof(double));
    hinterp = (double *)malloc((ndepths)*sizeof(double));
    for (i=0; i<ndepths; i++) { // Set y values to interpolate at
	yinterp[i] = c[i].y;
    }
    InterpolateArray(ya, fa, n, yinterp, finterp, ndepths, 0, 0.0, 1, 0.0);
    InterpolateArray(ya, ha, n, yinterp, hinterp, ndepths, 0, 0.0, 1, 0.0);
    for (i=0; i<ndepths; i++) { // Set u values in c based on interpolation
	c[i].f[k] = finterp[i];
	c[i].h[k] = hinterp[i];
    }
    free(ya);
    free(fa);
    free(ha);
    free(yinterp);
    free(finterp);
    free(hinterp);
}

