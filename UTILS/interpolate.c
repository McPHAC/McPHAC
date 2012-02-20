#include "includes.h"

void interpolate(double* xa, double* ya, int32_t n, double x, double* y, double* dy) {
	double *ypp, ypval, yppval;

	ypp = (double *) malloc((n)*sizeof(double));
	cspline_set(n, xa, ya, 0, 0.0, 0, 0.0, ypp);
	cspline_val(n, xa, ya, ypp, x, y, &ypval, &yppval);
	*dy = 0.0; // TODO: dy presently isn't supported, so perhaps remove arg
	free(ypp);
}
