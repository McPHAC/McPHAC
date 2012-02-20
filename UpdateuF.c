#include "includes.h"

void UpdateuF(int k, matrix_type M, column_type* c, int ndepths) {
	// Calculates P in the Feautrier solution to the radiative transfer

	double   *tmp1, *tmp2, *tmpu;
	int      i, j;

	tmp1 = dvector(NMU);
	tmp2 = dvector(NMU);
	tmpu = dvector(NMU);

	MatrixVectorMultiply(tmpu, c[ndepths-1].BFtmpinv, c[ndepths-1].QFtmp, NMU);
	for (j=0; j<NMU; j++) {
		c[ndepths-1].u[j][k]= tmpu[j];
	}
	for (i=(ndepths-2); i>=0; i--) {
		CalcCF(M.CF, i, k, c, ndepths);
		MatrixVectorMultiply(tmp1, M.CF, tmpu, NMU);
		for (j=0; j<NMU; j++) {
			tmp2[j] = c[i].QFtmp[j] - tmp1[j];
		}
		MatrixVectorMultiply(tmpu, c[i].BFtmpinv, tmp2, NMU);
		for (j=0; j<NMU; j++) {
			c[i].u[j][k]= tmpu[j];
		}
	}

	free(tmp1);
	free(tmp2);
	free(tmpu);
}
