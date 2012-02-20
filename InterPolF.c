#include "includes.h"

void InterPolF(column_type * cnu, int ndepthsnu, double *nu, int k, int knu, column_type * c, int ndepths)
// Fills the full column array c with F values based on interpolation of cnu
{
	int i,s,n;
	double  *ya, *Fa, *yinterp, *Finterp;

	for (s=0; (s<ndepths)?(c[s].y<=cnu[ndepthsnu-1].y):0; s++); // s is index of first depth of c beyond the range of cnu
	n = ndepthsnu + ndepths-s;
	ya = (double *)malloc((n)*sizeof(double));
	Fa = (double *)malloc((n)*sizeof(double));
	for (i=0; i<ndepthsnu; i++) {
		ya[i] = cnu[i].logy;
		Fa[i] = cnu[i].F[knu];
	}
	for (i=s; (i+1)<ndepths; i++) {
		ya[ndepthsnu+i-s] = c[i].logy;
		Fa[ndepthsnu+i-s] = 0.0;
	}
	ya[ndepthsnu+(ndepths-1)-s] = c[ndepths-1].logy;
	Fa[ndepthsnu+(ndepths-1)-s] = Fa[ndepthsnu+(ndepths-1)-s-1];
	yinterp = (double *)malloc((ndepths)*sizeof(double));
	Finterp = (double *)malloc((ndepths)*sizeof(double));
	for (i=0; i<ndepths; i++) { // set y values to interpolate at
		yinterp[i] = c[i].logy;
	}
	InterpolateArray(ya, Fa, n, yinterp, Finterp, ndepths, 0, 0.0, 1, 0.0);
	for (i=0; i<ndepths; i++) {
		c[i].F[k] = Finterp[i];
	}

	free(ya);
	free(Fa);
	free(yinterp);
	free(Finterp);
}
