#include "includes.h"


double *dvector(int n) {
    double *v;
    v= (double *)malloc((unsigned)n*sizeof(double));
    if (!v) {
        fprintf(stderr,"Failed to make Allocation of pointer to double...n=%d\n", n);
        exit(1);
    }
    memset(v, 0, n*sizeof(double)); /*zero the vector*/
    return(v);
}
