#include "includes.h"

column_type * AllocateColumns(int ndepths, int nmu, int nfreq) {
    // Allocates the column_type array

    column_type *c;
    int   i, j;


    c = (column_type *) malloc(ndepths*sizeof(column_type));
    for (i=0; i<ndepths; i++) {
        c[i].I=ddvector(nmu);
        c[i].u=ddvector(nmu);
        c[i].JanisT = ddvector(nmu);

        if (FEAUTRIER) {
	    c[i].BFtmp = ddvector(nmu);
	    c[i].BFtmpinv = ddvector(nmu);
	    c[i].QFtmp = dvector(nmu);
	}

        for (j=0; j<nmu; j++) {
            c[i].I[j]=dvector(nfreq);
            c[i].u[j]=dvector(nfreq);
            c[i].JanisT[j] = dvector(nfreq);

            if (FEAUTRIER) {
		c[i].BFtmp[j] = dvector(nmu);
		c[i].BFtmpinv[j] = dvector(nmu);
            }
        }

        c[i].tau = dvector(nfreq);
        c[i].dtau = dvector(nfreq);
        c[i].kappa = dvector(nfreq);
        c[i].k = dvector(nfreq);
        c[i].rho_opacity = dvector(nfreq);
        c[i].S = dvector(nfreq);
        c[i].J = dvector(nfreq);
        c[i].Jt = dvector(nfreq);
        c[i].h = dvector(nfreq);
        c[i].f = dvector(nfreq);
        c[i].q = dvector(nfreq);
        c[i].F = dvector(nfreq);
    }

    return(c);
}


