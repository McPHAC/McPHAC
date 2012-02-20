#include "includes.h"

void PrintOutOpacities(column_type *c) {
    // Prints the opacities at log y = -6, -4 and -2
    int i, k;
    int   y2, y4, y6;
    extern int NDepths;
    FILE  *fp;

    fp=fopen("OUT/opacities.dat", "w");
    fprintf(fp,"# E(kev)  log opacity @ logy=-2, -4, and -6\n");
    for (i=0; c[i].logy<-6 && i<NDepths; i++);
    y6=i;
    for (i=0; c[i].logy<-4 && i<NDepths; i++);
    y4=i;
    for (i=0; c[i].logy<-2 && i<NDepths; i++);
    y2=i;

    for (k=0; k<NFREQ; k++) {
        fprintf(fp,"%7.3g  %7.3g  %7.3g  %7.3g\n", nu[k]*HPLANCK/ERGSPERKEV, log10(c[y2].k[k]),
                log10(c[y4].k[k]), log10(c[y6].k[k]));
    }

    fclose(fp);
}



