/********************
 * The McGill Planar Hydrogen Atmosphere Code (McPHAC)
 ********************/

// Include file defining parameters, structs, etc.
#include "includes.h"

// Main program
int main( int argc, char *argv[] ) {
	// Re-declare external variables
	extern int NDepths, NGaunt;
	extern double **gauntff, *nu, *lognu, *dnu;
	extern matrix_type   M;


	// Declare variables
	double *mu2, *dmu2;
	int i, j, ionekev, iteration, lastndepths, stable, multiplier;
	int k, ndepthsnu, ndepths, nmu, nfreq, knu, nmunu, nfreqnu;
	int ndepthst, nfreqt, nmut;
	int ndepthsy, nfreqy, nmuy;
	column_type *c, *ct, *cy;
	column_type *cnu;
	void  legendre_com_(); // Fortran program, so declare it here
	double T;
	double dlognu, maxfreq, minfreq;
	FILE  *fp, *fpint, *fpJ, *fp8;
	char   filename[100];
	double  maxtemperaturechange;
	int dampflag;

	// Set variables given through input arguments
	if ( argc != 15+1 ) { // +1 is for executable name
		printf( "usage: %s Teff gsurface mincol maxcol ndepths maxfactor ndepthsnu maxfactornu nmu nfreq maxfractempchange maxiter anist maxcoltau tguessbdycond\n", argv[0] );
		return 0;
	} else {
		int i=1;
		TEFF = pow(10,atof(argv[i])); i++;
		GSURFACE = atof(argv[i]); i++;
		MINCOL = atof(argv[i]); i++;
		MAXCOL = atof(argv[i]); i++;
		NDEPTHS = atoi(argv[i]); i++;
		MAXFACTOR = atoi(argv[i]); i++;
		NDEPTHSNU = atoi(argv[i]); i++;
		MAXFACTORNU = atoi(argv[i]); i++;
		NMU = atoi(argv[i]); i++;
		NFREQ = atoi(argv[i]); i++;
		MAXFRACTEMPCHANGE = atof(argv[i]); i++;
		MAXITER = atoi(argv[i]); i++;
		ANIST = atoi(argv[i]); i++;
		MAXCOLTAU = atof(argv[i]); i++;
		TGUESSBDYCOND = atof(argv[i]); i++;
	}

	// Allocate arrays
	nu = (double *) malloc(NFREQ*sizeof(double));
	lognu = (double *) malloc(NFREQ*sizeof(double));
	dnu = (double *) malloc(NFREQ*sizeof(double));
	mu = (double *) malloc(NMU*sizeof(double));
	dmu = (double *) malloc(NMU*sizeof(double));
	dmudnu = (double *) malloc(NMU*NFREQ*sizeof(double));
	mu2 = (double *) malloc(2*NMU*sizeof(double));
	dmu2 = (double *) malloc(2*NMU*sizeof(double));

	// Initialize variables
	extern_maxtemperaturechange = 1.0;

	nmu = 2*NMU;
	legendre_com_(&nmu, mu2, dmu2); // Gauss-Legendre quadrature abscissas and weights on [-1,1]
	nmu = NMU;
	for (j=0; j<nmu; j++) { // Only use values for [0,1]
		mu[j] = mu2[NMU+j];
		dmu[j] = dmu2[NMU+j];
	}

	minfreq = (FREQKT) ? log10(MINFREQKT*KBOLTZMAN*HZPER1KEVPHOTON/ERGSPERKEV*TEFF) : MINFREQ;
	maxfreq = (FREQKT) ? log10(MAXFREQKT*KBOLTZMAN*HZPER1KEVPHOTON/ERGSPERKEV*TEFF) : MAXFREQ;
	dlognu = (maxfreq-minfreq)/(1.0*NFREQ); // Taking frequency values to be the centers of log spaced bins

	ionekev = 0;
	for (k=0; k<NFREQ; k++) { // Find k corresponding to 1 keV photon
		lognu[k] = minfreq + dlognu*k;
		nu[k] = pow(10.0, lognu[k]);
		if (nu[k]<HZPER1KEVPHOTON) {
			ionekev = k;
		}
		dnu[k] = nu[k]*dlognu*LNTEN;
	}

	ndepths = NDEPTHS*MAXFACTOR;
	nmu = (USEINTERPOL) ? 1 : NMU; // No angle dependence in the temperature correction procedure,
	// so only one angle point is needed if rad. transfer uses separate column array
	nfreq = NFREQ;
	c = GetColumnsLog(NULL, ndepths, nmu, nfreq, &ndepths, &nmu, &nfreq, 1);
	NDepths = ndepths; // In case the external variable NDepths still appears anywhere

	// Print out the initial structure and opacities at 1keV
	fp = fopen("OUT/InitialkR.dat", "w");
	fprintf(fp,"#logy, logT, rho, Pressure, kR, kappa (@ 1 keV), k (@ 1 keV)\n");
	for (i=0; i<NDepths; i++) {
		fprintf(fp,"%f  %f  %g  %g  %g %g  %g\n", c[i].logy, c[i].logT,
				c[i].rho, c[i].pressure, c[i].kR, c[i].kappa[ionekev],
				c[i].k[ionekev]);
	}
	fclose(fp);

	// Print out the Gaunt factor
	fp = fopen("OUT/GauntFactor.dat", "w");
	fprintf(fp,"#logy, nu, T, g2, u, gff\n");
	for (k=0; k<NFREQ; k++) {
		for (i=0; i<NDepths; i++) {
			T = pow(10.0,c[i].logT);
			fprintf(fp,"%f  %g  %g  %g  %g  %g\n", c[i].logy, nu[k], T,
					RY*ERGPEREV/T/KBOLTZMAN, nu[k]*HPLANCK/T/KBOLTZMAN,
					GauntFactor(nu[k], T));
		}

	}
	fclose(fp);

	// Proceeding to the calculation
	iteration = 1;
	stable = 0;
	dampflag = 0;
	maxtemperaturechange = 1.0;
	do {
		stable = (fabs(maxtemperaturechange)<MAXFRACTEMPCHANGE || iteration>MAXITER) ? 1 : 0;
		do {
			// File for outputting P
			sprintf(filename,"OUT/PymuInterp.%d.dat", iteration);
			fpint = fopen(filename, "w");
			fprintf(fpint, "#y mu nu tau P     P y mu nu tau\n");
			// File for outputting J calculated in the temperature correction
			sprintf(filename,"OUT/JT.%d.dat", iteration);
			fpJ = fopen(filename, "w");
			fprintf(fpJ, "#y mu nu J JT P tau\n");

			fprintf(stderr,"\n\n\nIteration %d\n\n", iteration);

			// Radiative Transfer Calculation
			nmunu = NMU;
			nfreqnu = nfreq;
			if (FEAUTRIER) {
				AllocateMatrix(&M, 1, nmunu);
			} else {
				ndepthsnu = NDEPTHSNU;
				AllocateMatrix(&M, ndepthsnu, 1);
			}
			for (multiplier=1; multiplier<=MAXFACTORNU; multiplier*=2) {
				ndepthsnu = (USEINTERPOL) ? NDEPTHSNU*multiplier : ndepths;
				fprintf(stderr, "Radiative transfer ndepthsnu=%d freq: ", ndepthsnu);

				// File for outputting P
				sprintf(filename,"OUT/Pymu.%d.%d.dat", ndepthsnu, iteration);
				fp = fopen(filename, "w");
				fprintf(fp, "#y mu nu tau P     P y mu nu tau\n");

				// File for outputting I at surface
				sprintf(filename,"OUT/EmergentSpectrum.%d.%d.dat", ndepthsnu, iteration);
				fp8 = fopen(filename, "w");
				fprintf(fp8, "#nu mu I\n");

				for (k=0; k<NFREQ; k++) {
					fprintf(stderr, "%d...", k);
					if (USEINTERPOL) {
						knu = 0;
						cnu = GetColumnsNu(c, ndepths, nmunu, k, 0, &ndepthsnu); // Generates cnu from c
					} else {
						knu = k;
						cnu = c;
					}

					if (FEAUTRIER) {
						CalcBF(cnu[0].BFtmp, 0, knu, cnu, ndepthsnu);
						CalcQF(cnu[0].QFtmp, 0 , k, knu, cnu, ndepthsnu);
						for (i=1; i<ndepthsnu; i++) {
							CalcAF(M.AF, i, knu, cnu, ndepthsnu);
							CalcBF(M.BF, i, knu, cnu, ndepthsnu);
							CalcCF(M.CF, i-1, knu, cnu, ndepthsnu);
							InvertMatrix(cnu[i-1].BFtmpinv, cnu[i-1].BFtmp, NMU);
							CalcBFtmp(cnu[i].BFtmp, M.AF, M.BF, cnu[i-1].BFtmpinv, M.CF);
							CalcQF(M.QF, i, k, knu, cnu, ndepthsnu);
							CalcQFtmp(cnu[i].QFtmp, M.AF, cnu[i-1].BFtmpinv, M.QF, cnu[i-1].QFtmp);
						}
						InvertMatrix(cnu[ndepthsnu-1].BFtmpinv, cnu[ndepthsnu-1].BFtmp, NMU);
						UpdateuF(knu, M, cnu, ndepthsnu);
						if (OUTPUT) {
							for (j=NMU; j<NMU; j++) {
								for (i=0; i<ndepthsnu; i++) {
									fprintf(fp, "%9.2e %9.2e %9.2e %9.2e %9.5e      %9.5e %9.2e %9.2e %9.2e %9.2e\n",
											cnu[i].y, mu[j], nu[knu], cnu[i].tau[knu], cnu[i].u[j][knu],
											cnu[i].u[NMU-1-j][knu], cnu[i].y, mu[NMU-1-j], nu[knu], cnu[i].tau[knu]);
								}
								fprintf(fp, "\n\n");
							}
						}
					} else {
						InitW(M.W, ndepthsnu); // Zero the W[ndepthsnu][ndepthsnu] matrix, adding -1 to the diagonal

						memset(M.Q, 0, ndepthsnu*sizeof(double));

						CalcU(M.U, knu, cnu, ndepthsnu);
						for (j=0; j<NMU; j++) {
							CalcK(M.K, j, k, knu, cnu, ndepthsnu);
							CalcT(M.T, j, knu, cnu, ndepthsnu);
							InvertTridiagonalMatrix(M.TInverse,M.T,ndepthsnu); // See note on ineff. of this in func. def.
							CalcV(M.V, j, k, cnu, ndepthsnu);
							UpdateW(M.W, M.V, M.TInverse, M.U, ndepthsnu);
							UpdateQ(M.Q, M.V, M.TInverse, M.K, ndepthsnu);
						}
						CalcJbar(M.W, M.Q, cnu, ndepthsnu);

						for (j=0; j<NMU; j++) {
							CalcK(M.K, j, k, knu, cnu, ndepthsnu);
							CalcT(M.T, j, knu, cnu, ndepthsnu);
							InvertTridiagonalMatrix(M.TInverse,M.T,ndepthsnu); // See note on ineff. of this in func. def.
							Updateu(j, knu, M.TInverse, M.K, M.U, cnu, ndepthsnu);
							if (OUTPUT) {
								for (i=0; i<ndepthsnu; i++) {
									fprintf(fp, "%9.2e %9.2e %9.2e %9.2e %9.5e      %9.5e %9.2e %9.2e %9.2e %9.2e\n",
											cnu[i].y, mu[j], nu[knu], cnu[i].tau[knu], cnu[i].u[j][knu],
											cnu[i].u[NMU-1-j][knu], cnu[i].y, mu[NMU-1-j], nu[knu], cnu[i].tau[knu]);
								}
								fprintf(fp, "\n\n");
							}
						}
					}


					if (USEINTERPOL) {
						CalcMeanIntensityJ(cnu, dmu, knu, ndepthsnu, nmunu);
						CalcEddingtonFactors(cnu, mu, dmu, knu, ndepthsnu, nmunu);
						CalcFluxes(cnu, mu, dmu, knu, ndepthsnu, nmunu);
						InterPolEdd(cnu, knu, ndepthsnu, c, k, ndepths);
						OutputEddingtonFactors(cnu, iteration, k, nu[k], knu, ndepthsnu);
						InterPolJ(cnu, ndepthsnu, nu, k, knu, c, ndepths);
						InterPolF(cnu, ndepthsnu, nu, k, knu, c, ndepths);
						OutputSpectrum(cnu, fp8, knu, nu, k, mu, nmunu);
					} else {
						CalcMeanIntensityJ(c, dmu, k, ndepths, nmu);
						CalcEddingtonFactors(c, mu, dmu, k, ndepths, nmu);
						OutputEddingtonFactors(c, iteration, k, nu[k], k, ndepths);
						CalcSpecificIntensities(c, mu, k, ndepths, nmu);
					}
					fprintf(fp, "\n\n");

					// Output interpolated u
					if (OUTPUT && !USEINTERPOL) {
						for (j=0; j<NMU; j++) {
							for (i=0; i<NDepths; i++) {
								fprintf(fpint, "%9.2e %9.2e %9.2e %9.2e %9.5e     %9.5e %9.2e %9.2e %9.2e %9.2e\n",
										c[i].y, mu[j], nu[k], c[i].tau[k], c[i].u[j][k],
										c[i].u[NMU-1-j][k], c[i].y, mu[NMU-1-j], nu[k], c[i].tau[k]);
							}
							fprintf(fpint, "\n\n");
						}
						fprintf(fpint, "\n\n");
					}

					if (USEINTERPOL) {
						FreeColumns(cnu, ndepthsnu, nmunu, 1);
					}
				}
				fprintf(stderr, "\n");

				fclose(fp);
				fclose(fp8);
			}

			if (FEAUTRIER) {
				FreeMatrix(&M, 1, nmunu);
			} else {
				FreeMatrix(&M, ndepthsnu, 1);
			}

			if (USEINTERPOL) {
				CalcTeff(c, nu, dnu, ndepths, nfreq);
				OutputFluxes(c, nu, dnu, ndepths, nfreq, iteration, 0);
			} else {
				CalcFluxesOld(c,iteration,ionekev);
				OutputSpectrumOld(c, iteration);
			}

			// Output anisotropic scattering JanisT comparison with J
			if (OUTPUT && !USEINTERPOL) {
				for (k=0; k<NFREQ; k++) {
					for (j=0; j<NMU; j++) {
						for (i=0; i<NDepths; i++) {
							fprintf(fpJ, "%9.2e %9.2e %9.2e %11.4e %11.4e %11.4e %9.2e\n",
									c[i].y, mu[j], nu[k], c[i].J[k], c[i].JanisT[j][k], c[i].u[j][k],
									c[i].tau[k]);
						}
						fprintf(fpJ, "\n\n");
					}
				}
			}


			// Temperature correction procedure
			nmuy = (USEINTERPOL) ? 1 : NMU;
			nfreqy = NFREQ;
			ndepthsy = NDEPTHS;
			if (MAXFACTOR>1) cy = GetColumnsLog(c, ndepths, nmu, nfreq, &ndepthsy, &nmuy, &nfreqy, 0); // Stores y vals of multiplier=1 column
			for (multiplier=1; multiplier<=MAXFACTOR; multiplier*=2) {
				nmut = (USEINTERPOL) ? 1 : NMU;
				nfreqt = NFREQ;
				ndepthst = NDEPTHS*multiplier;
				if (multiplier==MAXFACTOR) {
					ct = c;
				} else {
					ct = GetColumnsLog(c, ndepths, nmu, nfreq, &ndepthst, &nmut, &nfreqt, 0);
					for (k=0; k<nfreq; k++) {
						InterPolEdd(c, k, ndepths, ct, k, ndepthst);
					}
				}

				AllocateMatrix(&M, ndepthst, nmut);
				InitW(M.W, ndepthst); // Zero the W[ndepthst][ndepthst] matrix, adding -1 to diagonal

				memset(M.Q, 0, ndepthst*sizeof(double));

				fprintf(stderr,"Temp. corr. freq: ");
				for (k=0; k<NFREQ; k++) {
					CalcUt(M.U, k, ct, ndepthst);
					CalcKt(M.K, 0, k, ct, ndepthst);
					CalcTt(M.T, 0, k, ct, ndepthst);
					fprintf(stderr,"%d...", k);
					InvertTridiagonalMatrix(M.TInverse,M.T,ndepthst); // See note on ineff. of this in func. def.
					CalcVt(M.V, 0, k, ct, ndepthst);
					UpdateW(M.W, M.V, M.TInverse, M.U, ndepthst);
					UpdateQ(M.Q, M.V, M.TInverse, M.K, ndepthst);
				}
				fprintf(stderr,"\n");
				CalcJtbar(M.W, M.Q, ct, ndepthst);

				fprintf(stderr,"Temp. corr. freq: ");
				for (k=0; k<NFREQ; k++) {
					CalcUt(M.U, k, ct, ndepthst);
					CalcKt(M.K, 0, k, ct, ndepthst);
					CalcTt(M.T, 0, k, ct, ndepthst);
					fprintf(stderr,"%d...", k);
					InvertTridiagonalMatrix(M.TInverse,M.T,ndepthst); // See note on ineff. of this in func. def.
					UpdateJt(k, M.TInverse, M.K, M.U, ct, ndepthst, iteration, ionekev);
				}
				fprintf(stderr,"\n");
				FreeMatrix(&M, ndepthst, nmut);
				maxtemperaturechange = DeltaT(ct, ndepthst, iteration, maxtemperaturechange, &dampflag);
				fprintf(stderr,"Max Temperature Change: %f\n", maxtemperaturechange);
				PrintInterpolT(c, ndepths, ct, ndepthst, iteration);
				CalcFluxest(ct, nu, dnu, ndepthst, nfreqt); // uses F = 4pi*d/dtau(f_nu*J_nu)
				CalcTeff(ct, nu, dnu, ndepthst, nfreqt);
				OutputFluxes(ct, nu, dnu, ndepthst, nfreqt, iteration, 1);

				if (multiplier!=MAXFACTOR) {
					FreeColumns(ct, ndepthst, nmut, nfreqt);
				}
			}
			if (MAXFACTOR>1) FreeColumns(cy, ndepthsy, nmuy, nfreqy);

			iteration++;

			fclose(fpint);
			fclose(fpJ);
			if (fabs(maxtemperaturechange)>MAXFRACTEMPCHANGE) {
				fprintf(stderr,"Re-Adjusting Atmosphere\n");
				lastndepths = ndepths;
				c = GetColumnsLog(c, ndepths, nmu, nfreq, &ndepths, &nmu, &nfreq, 1);
			}
			extern_maxtemperaturechange = maxtemperaturechange;
		} while (fabs(maxtemperaturechange)>MAXFRACTEMPCHANGE && iteration<MAXITER);
	} while (stable == 0);

	for (i=0; i<NGaunt; i++) {
		free(gauntff[i]);
	}
	free(gauntff);

	// Free arrays
	free(nu);
	free(lognu);
	free(dnu);
	free(mu);
	free(dmu);
	free(dmudnu);
	free(mu2);
	free(dmu2);

	return 0;
}
