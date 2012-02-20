#include "includes.h"

void InterpolateArray(double* xdata, double* ydata, int ndata, double* xinterp, double* yinterp, int ninterp, int ibcbeg, double ybcbeg, int ibcend, double ybcend) {
	// Interpolates data at an array of values

	double *ypp, ypval, yppval, x, y;
	int i;

	ypp = (double *) malloc((ndata)*sizeof(double));
	cspline_set(ndata, xdata, ydata, ibcbeg, ybcbeg, ibcend, ybcbeg, ypp);

	for (i=0; i<ninterp; i++) {
		x = xinterp[i];
		cspline_val(ndata, xdata, ydata, ypp, x, &y, &ypval, &yppval);
		yinterp[i] = y;
	}

	free(ypp);
}
