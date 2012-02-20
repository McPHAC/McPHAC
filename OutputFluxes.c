#include "includes.h"


void OutputFluxes(column_type  *c, double *nu, double *dnu, int ndepths, int nfreq, int iteration, int tflag) {
	// Outputs various fluxes

	FILE  *fp, *fpF, *fpFnorm;
	char   filename[100];
	double  *nunorm, *lognunorm, *nunormint, *lognunormint, dlognunorm, *Fnorm, *logFnorm, *Fnormint, *logFnormint;
	int i, k, r;

	if (tflag) {
		sprintf(filename,"OUT/Fluxest.%d.%d.dat", ndepths, iteration);
	} else {
		sprintf(filename,"OUT/Fluxes.%d.%d.dat", ndepths, iteration);
	}
	fp=fopen(filename,"w");
	fprintf(fp,"#y   T   Teff\n");

	for (i=0; i<ndepths; i++) {
		fprintf(fp,"%e   %e   %e\n", c[i].y, pow(10.0,c[i].logT), c[i].Teff);
	}

	if (tflag) {
		sprintf(filename,"OUT/SurfaceFluxest.%d.%d.dat", ndepths, iteration);
	} else {
		sprintf(filename,"OUT/SurfaceFluxes.%d.%d.dat", ndepths, iteration);
	}
	fpF=fopen(filename,"w");
	fprintf(fpF,"#nu   flux\n");

	for (k=0; k<nfreq; k++) {
		fprintf(fpF,"%e   %e\n", nu[k], c[0].F[k]);
	}

	if (tflag) {
		sprintf(filename,"OUT/SurfaceFluxesNormt.%d.%d.dat", ndepths, iteration);
	} else {
		sprintf(filename,"OUT/SurfaceFluxesNorm.%d.%d.dat", ndepths, iteration);
	}
	fpFnorm=fopen(filename,"w");
	fprintf(fpFnorm,"#h*nu/k/Teff   F*nu/sigma/Teff^4   Teff   F_logint*nu/sigma/Teff^4\n");

	nunorm = (double *)malloc((nfreq)*sizeof(double));
	Fnorm = (double *)malloc((nfreq)*sizeof(double));
	lognunorm = (double *)malloc((nfreq)*sizeof(double));
	logFnorm = (double *)malloc((nfreq)*sizeof(double));
	for (k=0; k<nfreq; k++) {
		nunorm[k] = nu[k]*HPLANCK/KBOLTZMAN/c[0].Teff;
		Fnorm[k] = c[0].F[k]*nu[k]/SIGMA/pow(c[0].Teff, 4.0);
		lognunorm[k] = log10(nunorm[k]);
		logFnorm[k] = log10(Fnorm[k]);
	}

	nunormint = (double *)malloc((NFREQINTERP)*sizeof(double));
	Fnormint = (double *)malloc((NFREQINTERP)*sizeof(double));
	lognunormint = (double *)malloc((NFREQINTERP)*sizeof(double));
	logFnormint = (double *)malloc((NFREQINTERP)*sizeof(double));
	dlognunorm = (MAXFREQINTERP-MINFREQINTERP)/(1.0*NFREQINTERP);
	for (r=0; r<NFREQINTERP; r++) {
		lognunormint[r] = MINFREQINTERP + dlognunorm*r;
		nunormint[r] = pow(10.0, MINFREQINTERP + dlognunorm*r);
	}
	InterpolateArray(nunorm, Fnorm, nfreq, nunormint, Fnormint, NFREQINTERP, 0, 0.0, 0, 0.0);
	InterpolateArray(lognunorm, logFnorm, nfreq, lognunormint, logFnormint, NFREQINTERP, 0, 0.0, 0, 0.0);

	for (r=0; r<NFREQINTERP; r++) {
		nunormint[r] = pow(10.0, MINFREQINTERP + dlognunorm*r);
		if ((nunormint[r]<nu[0]*HPLANCK/KBOLTZMAN/c[0].Teff) ||
				(nunormint[r]>nu[nfreq-1]*HPLANCK/KBOLTZMAN/c[0].Teff)) {
			fprintf(fpFnorm,"%e   %e   %e   %e\n", nunormint[r], 1e-10, c[0].Teff, 1e-10);
		} else {
			fprintf(fpFnorm,"%e   %e   %e   %e\n", nunormint[r], Fnormint[r], c[0].Teff, pow(10.0,logFnormint[r]));
		}
	}


	fclose(fp);
	fclose(fpF);
	fclose(fpFnorm);
}


