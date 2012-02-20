#include "includes.h"


double DeltaT(column_type *c, int ndepths, int iteration, double prevmaxtemperaturechange, int *dampflag) {
	// Calculate the temperature correction at each depth from Jt

	char   filename[FILENAME];
	int  i, k;
	double maxtemperaturechange, wasT, nowT;
	double  sumnumerator, sumdenomenator,T;
	double  sumf, sumsqrf, *avgf, *sigmaf, f;
	extern double *nu, *dnu;
	FILE  *fp;
	int  maxtchangedepth, nf;

	double  test;
	double  Jnu, bbTnu;

	avgf = (double *) malloc (ndepths*sizeof(double));
	sigmaf = (double *) malloc (ndepths*sizeof(double));
	maxtemperaturechange = 0;

	for (i=0; i<ndepths; i++) {
		sumnumerator=0;
		sumf=sumsqrf=0;
		nf=0;
		test=sumdenomenator=0;
		T = pow(10.0,c[i].logT);
		for (k=0; k<NFREQ; k++) {
			sumdenomenator += c[i].kappa[k]*dBdT(nu[k],T)*dnu[k];
			Jnu = c[i].Jt[k];
			bbTnu = Bnu(nu[k], T);
			f= (Jnu-bbTnu)/bbTnu;
			sumf += f;
			sumsqrf += f*f;
			nf++;
			sumnumerator += c[i].kappa[k] * (Jnu - bbTnu) * dnu[k];
		}

		c[i].deltaT = sumnumerator/sumdenomenator;
		avgf[i]=sumf/nf;
		sigmaf[i] = sqrt(sumsqrf/nf - avgf[i]*avgf[i]);
	}

	wasT = 0.0;
	nowT = 0.0;
	maxtchangedepth = -1;
	sprintf(filename,"OUT/Temperatures.%d.%d.dat", ndepths, iteration);
	fprintf(stderr,"Writing temperatures to %s\n", filename);
	fp= fopen(filename, "w");
	fprintf(fp,"#d, logy, oldT newT deltaT deltaT/T, avgf, fabs(avgf/sigmaf)\n");
	for (i=0; i<ndepths; i++) {
		T = pow(10.0,c[i].logT);
		if (fabs(c[i].deltaT/T)>maxtemperaturechange) {
			maxtemperaturechange = fabs(c[i].deltaT/T);
			maxtchangedepth=i;
			wasT = c[i].logT;
			nowT = log10(T + c[i].deltaT);
		}
	}
	*dampflag = ((maxtemperaturechange/prevmaxtemperaturechange>MAXTCHANGERATIO) && (iteration>=MINITERDAMP)) ? 1 : *dampflag;
	for (i=0; i<ndepths; i++) {
		T = pow(10.0,c[i].logT);
		c[i].logT = (DAMPDT && dampflag) ?
				log10(T + DTFRAC*c[i].deltaT) : log10(T + c[i].deltaT);
		fprintf(fp,"%d  %f  %f  %f %f  %f     %g  %g\n", i, c[i].logy, T, pow(10.0,c[i].logT), c[i].deltaT, c[i].deltaT/T,
				avgf[i], fabs(avgf[i])/sigmaf[i]);
		c[i].deltaT=0;
	}
	fclose(fp);
	fprintf(stderr,"Output Temperature information: %s\n", filename);
	fprintf(stderr,"Maximum computed (not necessarily applied) fractional temperature change at depth %d (of %d): Was %f, now %f\n",
			maxtchangedepth, ndepths, wasT, nowT);
	free(avgf);
	free(sigmaf);
	return(maxtemperaturechange);
}
