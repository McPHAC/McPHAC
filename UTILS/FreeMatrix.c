#include "includes.h"

void FreeMatrix(matrix_type *L, int ndepths, int nmu) {
    // Frees the matrix_type

    int  i, j;

    for (i=0; i<ndepths; i++) {
        free(L->T[i]);
        free(L->TInverse[i]);
        free(L->U[i]);
        free(L->V[i]);
        free(L->W[i]);
    }
    free(L->K);
    free(L->Q);
    free(L->U);
    free(L->T);
    free(L->TInverse);
    free(L->V);
    free(L->W);

    for (j=0; j<nmu; j++) {
        free(L->AF[j]);
        free(L->BF[j]);
        free(L->CF[j]);
    }
    free(L->AF);
    free(L->BF);
    free(L->CF);
    free(L->QF);
}
