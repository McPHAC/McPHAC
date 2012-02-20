#include "includes.h"

void InvertTridiagonalMatrix(double** ans, double** A, int n) {
    // Inverts a tri-diagonal matrix

    void  dgtsv_();  // Fortran program, so declare it here
    int nn, nrhs, ldb, info;
    double *dl, *dd, *du, *b;

    int   i, j, k;

    // use a LAPACK linear equation solver to solve A*X=I so that X is Ainverse
    nn = n;
    nrhs = n;
    ldb = n;
    dl = (double *)malloc((nn-1)*sizeof(double));
    dd = (double *)malloc(nn*sizeof(double));
    du = (double *)malloc((nn-1+1)*sizeof(double));
    b = (double *)malloc((ldb*nn+1)*sizeof(double));
    dd[0] = A[0][0];
    du[0] = A[0][1];
    for (i=1; i<n; i++) {
	dl[i-1] = A[i][i-1];
	dd[i] = A[i][i];
	du[i] = A[i][i+1];
    }
    dl[nn-2] = A[nn-1][nn-2];
    dd[nn-1] = A[nn-1][nn-1];
    for (j=0; j<n; j++) {
	for (k=0; k<n; k++) {
	    b[j*n+k] = (j==k) ? 1.0 : 0.0;
	}
    }
    dgtsv_(&nn, &nrhs, dl, dd, du, b, &ldb, &info);
    for (j=0; j<n; j++) {
	for (i=0; i<n; i++) {
	    ans[i][j] = b[j*n+i];
	}
    }
    free(dl);
    free(dd);
    free(du); 
    free(b);
}

