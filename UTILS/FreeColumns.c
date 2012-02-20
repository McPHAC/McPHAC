#include "includes.h"

void
FreeColumns(column_type* c, int ndepths, int nmu, int nfreq) {
    int  i, j;

    for (i=0; i<ndepths; i++) {
        for (j=0; j<nmu; j++) {
            free(c[i].I[j]);
            free(c[i].u[j]);
            free(c[i].JanisT[j]);
	    if (FEAUTRIER) {
		free(c[i].BFtmp[j]);
		free(c[i].BFtmpinv[j]);
	    }
        }
        free(c[i].u);

        free(c[i].I);
        free(c[i].tau);
        free(c[i].dtau);
        free(c[i].kappa);
        free(c[i].k);
        free(c[i].S);
        free(c[i].J);
        free(c[i].JanisT);
	if (FEAUTRIER) {
	    free(c[i].BFtmp);
	    free(c[i].BFtmpinv);
            free(c[i].QFtmp);
        }
        free(c[i].Jt);
        free(c[i].h);
        free(c[i].f);
        free(c[i].q);
        free(c[i].F);
    }
}
