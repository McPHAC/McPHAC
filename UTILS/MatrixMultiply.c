#include "includes.h"

void MatrixMultiply(double** ans, double** A, double** B, int n) {
    // Multiplies two general square matrices

    void  dgemm_();  // Fortran program, so declare it here
    char transa, transb;
    int mm, nn, kk, lda, ldb, ldc;
    double alpha, beta, *a, *b, *c;
    int  i, j, k;
    double  sum;

    if (USEBLAS) {
        transa = 'n'; // do not transpose matrix a
        transb = 'n'; // do not transpose matrix b
        mm = n; // number of rows of matricies a and c
        nn = n; // number of columns of matrices b and c
        kk = n; // number of rows of matrix b and number of columns of matrix a
        alpha = 1.0; // scalar to multiply a*b by
        beta = 0.0; // scalar to multiply c by before adding alpha*a*b
        lda = mm;
        ldb = kk;
        ldc = mm;
        a = (double *)malloc(lda*kk*sizeof(double));
        b = (double *)malloc(ldb*nn*sizeof(double));
        c = (double *)malloc(ldc*nn*sizeof(double));
	for (k=0; k<kk; k++) {
	    for (i=0; i<mm; i++) {
	       a[k*mm+i] = A[i][k];
	    }
	}
  	for (j=0; j<nn; j++) {
	    for (k=0; k<kk; k++) {
	       b[j*kk+k] = B[k][j];
	    }
	}
        dgemm_(&transa,&transb,&mm,&nn,&kk,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
  	for (j=0; j<nn; j++) {
	    for (i=0; i<mm; i++) {
	       ans[i][j] = c[j*mm+i];
	    }
	}
        free(a);
        free(b);
        free(c);
    } else {
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		sum=0;
		for (k=0; k<n; k++) {
		    sum += A[i][k] * B[k][j];
		}
		ans[i][j]=sum;
	    }
	}
    }
}


