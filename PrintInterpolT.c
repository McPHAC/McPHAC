#include "includes.h"

void PrintInterpolT(column_type *c, int ndepths, column_type *ct, int ndepthst, int iteration)
// Prints temperature profile stored in ct interpolated to column depths of c
{
    char   filename[FILENAME];
    FILE  *fp;
    double *logT, *logTt, *logy, *logyt;
    double *T, *Tt, *y, *yt;
    int i;

    logT = (double *)malloc((ndepths)*sizeof(double));
    logy = (double *)malloc((ndepths)*sizeof(double));
    logTt = (double *)malloc((ndepthst)*sizeof(double));
    logyt = (double *)malloc((ndepthst)*sizeof(double));
    T = (double *)malloc((ndepths)*sizeof(double));
    y = (double *)malloc((ndepths)*sizeof(double));
    Tt = (double *)malloc((ndepthst)*sizeof(double));
    yt = (double *)malloc((ndepthst)*sizeof(double));

    for (i=0; i<ndepths; i++) {
        logy[i] = log10(c[i].y);
        y[i] = c[i].y;
    }

    for (i=0; i<ndepthst; i++) {
        logyt[i] = log10(ct[i].y);
        logTt[i] = ct[i].logT;
        yt[i] = ct[i].y;
        Tt[i] = pow(10.0,ct[i].logT);
    }

    InterpolateArray(logyt, logTt, ndepthst, logy, logT, ndepths, 0, 0.0, 0, 0.0);
    InterpolateArray(yt, Tt, ndepthst, y, T, ndepths, 0, 0.0, 0, 0.0);

    sprintf(filename,"OUT/TempProfile.%d.%d.dat", ndepthst, iteration);
    fprintf(stderr,"Writing interpolated temperatures to %s\n\n", filename);
    fp= fopen(filename, "w");
    fprintf(fp,"#y c.T ct.T ct.Tlogint\n");

    for (i=0; i<ndepths; i++) {
        fprintf(fp,"%.15e  %.15e  %.15e  %.15e\n", c[i].y, pow(10.0,c[i].logT), T[i], pow(10.0,logT[i]));
    }

    fclose(fp);
}
