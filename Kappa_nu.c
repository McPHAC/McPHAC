#include "includes.h"

double Kappa_nu(double rho, double T, double localnu) {
	// Returns free free opacity, in cm^2/g
	double ans, prefactor;

	prefactor = 4.0 * pow(ECHARGE,6.0)/ (3.0 * MELECTRON*HPLANCK* CCC) *
			sqrt(2.0*M_PI/(3.0*KBOLTZMAN*MELECTRON))/SQR(MPROTON  + MELECTRON);
	ans = (FREEFREE) ? prefactor * 1./sqrt(T) * rho * pow(localnu, -3) *
			(1.0 - exp(-HPOVERKB*localnu/(T))) * GauntFactor(localnu, T) : 0.0;
	if (GREY) ans=SIGMA_T/(MPROTON + MELECTRON); // Use a frequency indep. opacity if GREY flag set
	return(ans);
}
