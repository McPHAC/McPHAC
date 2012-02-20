#include "includes.h"

void OutputEddingtonFactors(column_type  *c, int iteration, int k, double nu, int knu, int ndepths) {
	// Outputs the Eddington factors

	int i;
	FILE  *fp;
	char   filename[FILENAME];

	if (k==0) {
		sprintf(filename,"OUT/EddingtonFactors.%d.%d.dat", ndepths, iteration);
		fp= fopen(filename, "w");
		fprintf(fp,"#f, h, f*J, h*J q y tau\n");
	} else {
		sprintf(filename,"OUT/EddingtonFactors.%d.%d.dat", ndepths, iteration);
		fp= fopen(filename, "a");
	}
	for (i=0; i<ndepths; i++) {
		fprintf(fp,"%e  %e  %e  %e  %e  %e  %e  %e\n", c[i].f[knu], c[i].h[knu], c[i].f[knu]*c[i].J[knu],
				c[i].h[knu]*c[i].J[knu], c[i].q[knu], c[i].y, c[i].tau[knu], nu);
	}
	fprintf(fp,"\n\n");
	fclose(fp);
}
