#include "includes.h"

void AllocateMatrix(matrix_type *L, int ndepths, int nmu) {
	// Allocates the matrix_type

	int  i, j;

	L->K=dvector(ndepths);
	L->Q=dvector(ndepths);
	L->U=ddvector(ndepths);
	L->T=ddvector(ndepths);
	L->TInverse=ddvector(ndepths);
	L->V=ddvector(ndepths);
	L->W=ddvector(ndepths);

	for (i=0; i<ndepths; i++) {
		L->T[i]=dvector(ndepths);
		L->TInverse[i] = dvector(ndepths);
		L->U[i]=dvector(ndepths);
		L->V[i]=dvector(ndepths);
		L->W[i]=dvector(ndepths);
	}

	L->AF=ddvector(nmu);
	L->BF=ddvector(nmu);
	L->CF=ddvector(nmu);
	L->QF=dvector(nmu);

	for (j=0; j<nmu; j++) {
		L->AF[j]=dvector(nmu);
		L->BF[j]=dvector(nmu);
		L->CF[j]=dvector(nmu);
	}

	fprintf(stderr,"Done Allocating Memory to Matrices KUTVW and TInverse....\n");
}

