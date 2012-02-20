#include "includes.h"

void Updateu(int j, int k, double** TInverse, double* K, double **U, column_type* c, int ndepths) {
	// Calculates P in the Rybicki method in the radiative transfer

	double   *tmp1, **tmp2, *tmp3, u;
	int      i, s;

	tmp1 = dvector(ndepths);
	tmp2 = ddvector(ndepths);
	tmp3 = dvector(ndepths);
	for (i=0; i<ndepths; i++) {
		tmp2[i]=dvector(ndepths);
	}

	MatrixVectorMultiply(tmp1, TInverse, K, ndepths);
	MatrixDiagonalMatrixMultiply(tmp2,TInverse,U, ndepths);


	for (i=0; i<ndepths; i++) {
		tmp3[i]=0.0;
		for (s=0; s<ndepths; s++) {
			tmp3[i] += tmp2[i][s] * c[s].Jbar;
		}
		u=tmp1[i] - tmp3[i];
		c[i].u[j][k]= u;
	}

	for (i=0; i<ndepths; i++) {
		free(tmp2[i]);
	}

	free(tmp1);
	free(tmp2);
	free(tmp3);
}
