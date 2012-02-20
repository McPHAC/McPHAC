#include "includes.h"

void OutputSpectrumOld(column_type *c, int iteration)
// Outputs the surface specific intensity as a function of angle and frequency
{
    FILE *fp;
    int j, k;
    char   filename[30];

    sprintf(filename,"OUT/EmergentSpectrum.%d.dat", iteration);
    fp = fopen(filename, "w");

    fprintf(fp, "#    nu\\mu ");
    for (j=0; j<NMU; j++){ // across all angles
	fprintf(fp, " %9.2e", mu[j]);
    }
    fprintf(fp, "\n");

    for (k=0; k<NFREQ; k++){ // for all frequencies considered
	fprintf(fp, "%9.2e ", nu[k]);
	for (j=0; j<NMU; j++){ // across all angles
            fprintf(fp, " %9.2e", c[0].I[j][k]);
	}
	fprintf(fp, "\n");
    }
    fclose(fp);
}
