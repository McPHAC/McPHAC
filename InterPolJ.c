#include "includes.h"

void InterPolJ(column_type * cnu, int ndepthsnu, double *nu, int k, int knu, column_type * c, int ndepths)
// Fills the full column array c with J values based on BB function and interpolation of cnu
{
	int i,s,n;
	double T;
	double  *ya, *Ja, *yinterp, *Jinterp;

	for (s=0; (s<ndepths)?(c[s].y<=cnu[ndepthsnu-1].y):0; s++); // s is index of first depth of c beyond the range of cnu
	n = ndepthsnu + ndepths-s;
	ya = (double *)malloc((n)*sizeof(double));
	Ja = (double *)malloc((n)*sizeof(double));
	for (i=0; i<ndepthsnu; i++) { // Set y, f, and h values in cnu
		ya[i] = cnu[i].y;
		Ja[i] = cnu[i].J[knu];
	}
	for (i=s; i<ndepths; i++) { // Set y, f, and h values beyond depth s
		ya[ndepthsnu+i-s] = c[i].y;
		T = pow(10.0, c[i].logT);
		Ja[ndepthsnu+i-s] = Bnu(nu[k],T);
	}
	yinterp = (double *)malloc((ndepths)*sizeof(double));
	Jinterp = (double *)malloc((ndepths)*sizeof(double));
	for (i=0; i<ndepths; i++) { // Set y values to interpolate at
		yinterp[i] = c[i].y;
	}
	InterpolateArray(ya, Ja, n, yinterp, Jinterp, ndepths, 0, 0.0, 0, 0.0);
	for (i=0; i<ndepths; i++) { // Set u values in c based on interpolation
		c[i].J[k] = Jinterp[i];
	}

	free(ya);
	free(Ja);
	free(yinterp);
	free(Jinterp);
}

