#include "includes.h"

column_type* GetColumnsLog(column_type *c, int ndepths, int nmu, int nfreq, int *newndepths, int *newnmu, int *newnfreq, int freeflag) {
    // Creates or updates the column struct array

    int  i, k, lastcol;
    extern int NDepths;
    column_type *newc;
    double  dlogy;
    double  maxtau;
    FILE  *fp;
    double  *logy, *newlogy, *logT, *newlogT, Tnorm;

    if (c==NULL) { // Doesn't exist, so create from scratch
        fp = fopen("OUT/Thermo.dat","w");
        if (fp==NULL) {
            fprintf(stderr,"Failed to open file OUT/Thermo.dat.  Does OUT directory exist?\nExiting...\n");
            exit(-1);
        }

        c = AllocateColumns(ndepths, nmu, nfreq);
        dlogy = 1.0*(MAXCOL-MINCOL)/(1.0*ndepths);
        maxtau = 0.0;
        while (maxtau<MAXCOLTAU) { // Ensure that many enough optical depths are considered at the highest frequency
            if (maxtau>0.0) {
                fprintf(stderr, "Largest optical depth at max frequency is %e ...increasing maxcol (log) by 10%%.\n", maxtau);
                maxtau = 0.0;
                dlogy += 0.1*fabs(MAXCOL)/(1.0*ndepths);
            }
            lastcol=0;
	    // Calculate initial guess at atmospheric structure
            for (i=0; i<ndepths; i++) {
                c[i].logy = MINCOL + dlogy * i;
                c[i].y = pow(10.0,c[i].logy);
                c[i].pressure = GSURFACE * c[i].y;
                if (i==0) {
                    c[i].logT = log10(TGUESSBDYCOND*TEFF);
                    c[ndepths-1].dy = 0.0; // There is no column after the deepest depth
                } else {
                    c[i-1].dy = c[i].y-c[i-1].y; // dy is the amount of column between the present and next depth
                    Tnorm = pow(10.0, c[i-1].logT)/TEFF;
                    c[i].logT = log10(TEFF*(Tnorm+3.0/16.0*c[i-1].kR*pow(Tnorm, -3.0)*c[i-1].dy));
                }
                c[i].rho = CalcRho(c[i].pressure, pow(10.0, c[i].logT));
                CalcOpacities(&c[i], 0, NFREQ, 0);
                CalckR(&c[i]);
                if (i==0) {
                    for (k=0; k<nfreq; k++) {
                        c[0].tau[k] = c[0].k[k]*c[0].y;
                        c[ndepths-1].dtau[k] = 0.0;
                    }
                } else {
                    for (k=0; k<nfreq; k++) {
                        c[i-1].dtau[k] = 0.5*(c[i].k[k]+c[i-1].k[k])*c[i-1].dy;
                        c[i].tau[k] = c[i-1].tau[k] + c[i-1].dtau[k];
                    }
                }
                if (i>0) maxtau += 0.5*(c[i].k[nfreq-1] + c[i-1].k[nfreq-1])*c[i-1].dy;
                if (maxtau<MAXCOLTAU) lastcol++;
            }
        }

        *newndepths = ndepths;
        *newnmu = nmu;
        *newnfreq = nfreq;

        fprintf(stderr,"New # of column depths: %d\n", *newndepths);

        fprintf(fp,"# logy logT pressure rho kR");
        for (i=0; i<*newndepths; i++) {
            fprintf(fp,"\n%f %f %g %g %g", c[i].logy, c[i].logT, c[i].pressure, c[i].rho, c[i].kR);
        }
        fprintf(fp,"\n");
        fclose(fp);
        NDepths = *newndepths;
        PrintOutTau(c);
        PrintOutOpacities(c);
        return (c);
    } else { // Exists, so update based on new temperature profile
        fprintf(stderr,"Resetting column depths..\n");

        if (freeflag) NDepths = *newndepths;
        newc = AllocateColumns(*newndepths, *newnmu, *newnfreq);

        logy = (double *) malloc((ndepths)*sizeof(double));
        logT = (double *) malloc((ndepths)*sizeof(double));
        for (i=0; i<ndepths; i++) {
            logy[i] = c[i].logy;
            logT[i] = c[i].logT;
        }
        newlogy = (double *) malloc((*newndepths)*sizeof(double));
        newlogT = (double *) malloc((*newndepths)*sizeof(double));
        dlogy = 1.0*(c[ndepths-1].logy-c[0].logy)/(1.0*(*newndepths-1));
        for (i=0; i<*newndepths; i++) {
            newlogy[i] = c[0].logy + dlogy*i;
        }
        InterpolateArray(logy, logT, ndepths, newlogy, newlogT, *newndepths, 0, 0.0, 0, 0.0);

        for (i=0; i<*newndepths; i++) {
            newc[i].y=pow(10.0, newlogy[i]);
            newc[i].logy = newlogy[i];
            newc[i].logT = newlogT[i];
            newc[i].pressure = GSURFACE * newc[i].y;
            newc[i].rho = CalcRho(newc[i].pressure, pow(10.0, newc[i].logT));
            CalcOpacities(&newc[i], 0, NFREQ, 0);
            CalckR(&newc[i]);
            if (i+1<*newndepths) {
                newc[i].dy = pow(10.0, newlogy[i+1])-pow(10.0, newlogy[i]);
            } else {
                newc[i].dy=0;
            }
        }
        for (k=0; k<nfreq; k++) {
            newc[0].tau[k] = newc[0].k[k]*newc[0].y;
        }
        for (i=0; i<*newndepths; i++) {
            if (i+1<*newndepths) {
                for (k=0; k<nfreq; k++) {
                    newc[i].dtau[k] = 0.5*(newc[i].k[k]+newc[i+1].k[k])*newc[i].dy;
                    newc[i+1].tau[k] = newc[i].tau[k] + newc[i].dtau[k];
                }
            } else {
                for (k=0; k<nfreq; k++) {
                    newc[i].dtau[k] = 0.0;
                }
            }
        }

        free(logy);
        free(logT);
        free(newlogy);
        free(newlogT);
        if (freeflag) FreeColumns(c, ndepths, nmu, nfreq);

        fprintf(stderr,"New # of column depths: %d\n", *newndepths);
        return (newc);
    }
}
