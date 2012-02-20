#include "includes.h"

void OutputSpectrum(column_type *c, FILE *fp, int knu, double *nu, int k, double *mu, int nmu)
// Outputs the surface specific intensity at surface as a function of angle for a given frequency
{
	int j;

	for (j=0; j<nmu; j++){
		fprintf(fp, "%9.4e %9.2e %11.4e\n", nu[k], mu[j], 2.0*c[0].u[j][knu]);
	}
	fprintf(fp, "\n\n");
}


