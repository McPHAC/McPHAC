#include "includes.h"

void
CalcBFtmp(double **BFtmp, double **AF, double **BF, double **prevBFtmpinv, double **CF) {
	// Calculates an intermediate matrix in the Feautrier solution to the radiative transfer
	double   **tmp1, **tmp2;
	int      j, jj;

	tmp1 = ddvector(NMU);
	tmp2 = ddvector(NMU);
	for (j=0; j<NMU; j++) {
		tmp1[j]=dvector(NMU);
		tmp2[j]=dvector(NMU);
	}

	MatrixMultiply(tmp1, prevBFtmpinv, CF, NMU);
	MatrixMultiply(tmp2, AF, tmp1, NMU);
	for (j=0; j<NMU; j++) {
		for (jj=0; jj<NMU; jj++) {
			BFtmp[j][jj] = BF[j][jj] - tmp2[j][jj];
		}
	}

	for (j=0; j<NMU; j++) {
		free(tmp1[j]);
		free(tmp2[j]);
	}
	free(tmp1);
	free(tmp2);
}
