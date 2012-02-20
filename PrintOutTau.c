#include "includes.h"

void PrintOutTau(column_type *c) {
	// Prints the column density at which the optical depth is 1, 2, 5, and 10, for each freq.

	extern int NDepths;
	extern  double  *nu;
	int  i, k;
	int  found1, found2, found5, found10;
	double   opacity;
	FILE   *fp;

	fp=fopen("OUT/tau.dat", "w");
	fprintf(fp,"#E(kev)  y(tau=1)  y(tau=2)  y(tau=5)  y(tau=10)");
	for (k=0; k<NFREQ; k++) {
		opacity=0;
		fprintf(fp,"\n%g ", nu[k]*HPLANCK/ERGSPERKEV);
		found1=found2=found5=found10=0;
		for (i=1; i<NDepths; i++) {
			opacity += 0.5*(c[i].k[k] +c[i-1].k[k])*c[i-1].dy;
			if (opacity>1 && !found1) {
				found1=1;
				fprintf(fp, "%g ", c[i].y);
			}
			if (opacity>2 && !found2) {
				found2=2;
				fprintf(fp, "%g ", c[i].y);
			}
			if (opacity>5 && !found5) {
				found5=2;
				fprintf(fp, "%g ", c[i].y);
			}
			if (opacity>10 && !found10) {
				found10=2;
				fprintf(fp, "%g ", c[i].y);
				break;
			}
		}
	}
	fclose(fp);
}

