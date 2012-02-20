#include "includes.h"

void InvertMatrix(double** ans, double** A, int n) {
	// Inverts a general matrix

	// TODO: Note that forming an explicit inverse here should be less efficient than
	// doing a solve wherever needed, so it might be worthwhile to change this

	void  dgesv_();  // Fortran program, so declare it here
	int nn, nrhs, lda, ldb, *ipiv, info;
	double *a, *b;

	int   i, j, k;

	// Use a LAPACK linear equation solver to solve A*X=I so that X is Ainverse
	nn = n;
	nrhs = n;
	lda = n;
	ldb = n;
	a = (double *)malloc(lda*nn*sizeof(double));
	b = (double *)malloc(ldb*nn*sizeof(double));
	ipiv = (int *)malloc(nn*sizeof(int));
	for (k=0; k<n; k++) {
		for (i=0; i<n; i++) {
			a[k*n+i] = A[i][k];
		}
	}
	for (j=0; j<n; j++) {
		for (k=0; k<n; k++) {
			b[j*n+k] = (j==k) ? 1.0 : 0.0;
		}
	}
	dgesv_(&nn, &nrhs, a, &lda, ipiv, b, &ldb, &info);
	for (j=0; j<n; j++) {
		for (i=0; i<n; i++) {
			ans[i][j] = b[j*n+i];
		}
	}
	free(a);
	free(b);
	free(ipiv);
}


