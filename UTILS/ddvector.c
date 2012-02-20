#include "includes.h"


double **ddvector(int n) {
    double **v;
    v= (double **)malloc((unsigned)n*sizeof(double *));
    if (!v) {
        fprintf(stderr,"Failed to make Allocation of pointer to pointer to double...n=%d\n", n);
        exit(1);
    }
    return(v);
}
