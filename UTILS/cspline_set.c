#include "includes.h"

void cspline_set(int32_t n, double* xa, double* ya, int32_t ibcbeg, double ybcbeg, int32_t ibcend, double ybcend, double* ypp) {
    void  spline_cubic_set_(); 
    spline_cubic_set_(&n, xa, ya, &ibcbeg, &ybcbeg, &ibcend, &ybcend, ypp);
}
