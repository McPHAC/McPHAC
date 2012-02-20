#include "includes.h"

void UpdateJt(int k, double** TInverse, double* K, double **U, column_type* c, int ndepths, int iteration, int ionekev) {
	// Calculates Jt in the temperature correction procedure

	double   *tmp1, **tmp2, *tmp3, u;
	int      i, s;
	FILE  *fp;
	char   filename[FILENAME];

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
			tmp3[i] += tmp2[i][s] * c[s].Jtbar;
		}
		u=tmp1[i] - tmp3[i];
		c[i].Jt[k]= u;
	}

	for (i=0; i<ndepths; i++) {
		free(tmp2[i]);
	}

	free(tmp1);
	free(tmp2);
	free(tmp3);

	sprintf(filename,"OUT/JvsJt.%d.dat", iteration);
	fp= fopen(filename, "w");
	fprintf(fp,"#J, Jt, Jt/J\n");
	for (i=0; i<ndepths; i++) {
		fprintf(fp,"%e  %e  %e\n", c[i].J[ionekev], c[i].Jt[ionekev], c[i].Jt[ionekev]/c[i].J[ionekev]);
	}
	fclose(fp);
}

