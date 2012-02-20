#include "includes.h"

double Sigma(double fioniz) {
    // Returns the scattering opacity (cm^2/g)
    double ans;
    ans = (THOMSON) ? SIGMA_T/(MPROTON + MELECTRON) : 0.0;

    return(ans);
}
