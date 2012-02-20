#include "includes.h"

void cspline_val(int32_t n, double* xa, double* ya, double* ypp, double x, double* y, double* ypval, double* yppval) {
    void  spline_cubic_val_();  /*this is a fortran program, so we'll declare it here*/
    spline_cubic_val_(&n, xa, ya, ypp, &x, y, ypval, yppval);
}

