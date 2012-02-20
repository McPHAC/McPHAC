#include "includes.h"

void CalcOpacities(column_type  * c, int kmin, int nk, int kk) {
	// Assigns opacities to a column type struct

	int  k;
	extern  double  *nu;
	double T, rho;

	rho=c->rho;
	T = pow(10.0,c->logT);
	c->sigma = Sigma(1.0); // Scattering cross section, cm^2/g

	for (k=kmin; k<kmin+nk; k++) {
		kk = (nk==1) ? kk : k; // input kk, if nk=1, specifies the k index of nu[k]
		c->kappa[k] = Kappa_nu(rho, T, nu[kk]); // True absorption opacity, cm^2/g
		c->k[k] = c->kappa[k] + c->sigma; // Total absorption opacity, cm^2/g
		c->rho_opacity[k] = c->sigma/c->k[k]; // Dimensionless scattering albedo
	}
}
