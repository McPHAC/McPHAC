#include "includes.h"

#define NINTERPOLATE  3  /*we use 2*NINTERPOLATE values for a polynomial interpretation*/
void CalcFluxesOld(column_type  *c, int  iteration, int ionekev) {
	// Old code for calculating and outputting various quantities related to energy flux,
	// used only when the INTERPOLATE flag is not set.

	int   i, j, k, r, kstart, kend, kk;
	double  flux, fluxanisT, totalFlux;
	extern double *mu, *dmu, *nu;
	FILE  *fp, *fpI, *fpE, *fpF, *fpFT, *fpFnorm;
	char   filename[30];
	double  x[NINTERPOLATE*2 +1], y[NINTERPOLATE*2+1], nunorm, dlognunorm, Fnorm, dFnorm;

	sprintf(filename,"OUT/Fluxes.%d.dat", iteration);
	fp=fopen(filename,"w");
	fprintf(fp,"# log Y  Flux,  log thermodynamic T (K), log Teff\n");

	sprintf(filename,"OUT/SpecificI.%d.dat", iteration);
	fpI = fopen(filename, "w");
	fprintf(fpI, "# y mu nu I IanisT tau\n");
	if (OUTPUT) {
		for (k=0; k<NFREQ; k++) {
			for (j=0; j<NMU; j++) {
				for (i=0; i<NDepths; i++) {
					fprintf(fpI, "%9.2e %9.2e %9.2e %11.4e %11.4e %9.2e\n",
							c[i].y, mu[j], nu[k], c[i].I[j][k], 0.0, c[i].tau[k]);
				}
				fprintf(fpI, "\n\n");
			}
		}
	}

	sprintf(filename,"OUT/SurfaceFluxes.%d.dat", iteration);
	fpF=fopen(filename,"w");
	fprintf(fpF,"nu   Flux\n");

	sprintf(filename,"OUT/SurfaceFluxesT.%d.dat", iteration);
	fpFT=fopen(filename,"w");
	fprintf(fpFT,"nu   Flux\n");

	sprintf(filename,"OUT/SurfaceFluxesNorm.%d.dat", iteration);
	fpFnorm=fopen(filename,"w");
	fprintf(fpFnorm,"h*nu/k/Teff   F*nu/sigma/Teff^4   Teff\n");

	sprintf(filename,"OUT/CalculationErrors.%d.dat", iteration);
	fpE=fopen(filename,"w");
	if (MESSAGES) fprintf(stderr,"Calculating Fluxes ... %s\n", filename);
	for (i=0; i<NDepths; i++) {
		totalFlux=0;
		for (k=0; k<NFREQ; k++) {
			flux=0;
			fluxanisT=0;
			for (j=0; j<NMU; j++) {
				flux += 4.0*M_PI*c[i].I[j][k] * mu[j]*dmu[j];
			}
			// Do a quick I+ > I- check
			for (j=NMU; j<NMU; j++) {
				if (c[i].I[j][k]< c[i].I[NMU-1-j][k]) {
					fprintf(fpE,"CALCULATION ERROR: I(-mu)[depth=%9.2e (%3d/%d)][mu=%f (%3d)][nu=%3d (%8.3g)](%g)>I(+mu)(%g).  Computational 	error.\n", c[i].y, i, NDepths-1, mu[j], j, k, nu[k], c[i].I[NMU-1-j][k], c[i].I[j][k]);
					if (c[i].I[j][k]/c[i].I[NMU-1-j][k] < 0.95*fabs(mu[j])) { // Prints to stderr if large error
						fprintf(stderr,"CALCULATION ERROR: I(-mu)[depth=%9.2e (%3d/%d)][mu=%f (%3d)][nu=%3d (%8.3g)](%g)>I(+mu)(%g).  Computational error.\n", c[i].y, i, NDepths-1, mu[j], j, k, nu[k], c[i].I[NMU-1-j][k], c[i].I[j][k]);
					}
				}
			}
			c[i].F[k]=flux;
			totalFlux += flux*dnu[k];
			if (i==0) {
				fprintf(fpF,"%e   %e\n", nu[k], flux);
				fprintf(fpFT,"%e   %e\n", nu[k], fluxanisT);
			}
		}
		c[i].totalFlux=totalFlux;
		c[i].Teff = pow(totalFlux/SIGMA,0.25);
		fprintf(fp,"%f   %g   %g   %g\n", c[i].logy, c[i].totalFlux, c[i].logT, log10(c[i].Teff));
	}

	for (r=0; r<NFREQINTERP; r++) {
		dlognunorm = (MAXFREQINTERP-MINFREQINTERP)/(1.0*NFREQINTERP);
		nunorm = pow(10.0, MINFREQINTERP + dlognunorm*r);
		if ((nunorm<nu[0]*HPLANCK/KBOLTZMAN/c[0].Teff) || (nunorm>nu[NFREQ-1]*HPLANCK/KBOLTZMAN/c[0].Teff)) {
			fprintf(fpFnorm,"%e   %e   %e\n", nunorm, 1e-10, c[0].Teff);
		} else {
			// Polynomial interpolation of normalized spectrum
			kk = 0;
			while (nu[kk]*HPLANCK/KBOLTZMAN/c[0].Teff<nunorm && kk<NFREQ-1) kk++;
			kstart = (kk-NINTERPOLATE>-1) ? kk-NINTERPOLATE : 0;
			kend = (kk+NINTERPOLATE<NFREQ) ? kk+NINTERPOLATE : NFREQ;
			if (kstart==0) {
				kend = 2*NINTERPOLATE;
			}
			if (kend==NFREQ) {
				kstart = NFREQ-2*NINTERPOLATE;
			}
			for (k=kstart; k<kend; k++) {
				x[k-kstart] = nu[k]*HPLANCK/KBOLTZMAN/c[0].Teff;
				y[k-kstart] = c[0].F[k]*nu[k]/SIGMA/pow(c[0].Teff, 4.0);
			}
			interpolate(x, y, kend-kstart, nunorm, &Fnorm, &dFnorm);
			fprintf(fpFnorm,"%e   %e   %e\n", nunorm, Fnorm, c[0].Teff);
		}
	}

	fseek ( fpE , 0 , SEEK_END );
	if ( ftell(fpE)==0) { // Remove error file if empty
		sprintf(filename,"OUT/CalculationErrors.%d.dat", iteration);
		remove(filename);
	} else fclose(fpE);
	fclose(fpF);
	fclose(fpFT);
	fclose(fpFnorm);
	fclose(fp);
	fclose(fpI);
	return;
}

#undef NINTERPOLATE

