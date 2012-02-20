#include "includes.h"

void
CalcQFtmp(double *QFtmp, double **AF, double **prevBFtmpinv, double *QF, double *prevQFtmp) {
    // Calculates a temporary vector used in the Feautrier solution to the radiative transfer
    double   *tmp1, *tmp2;
    int      j;

    tmp1 = dvector(NMU);
    tmp2 = dvector(NMU);

    MatrixVectorMultiply(tmp1, prevBFtmpinv, prevQFtmp, NMU);
    MatrixVectorMultiply(tmp2, AF, tmp1, NMU);
    for (j=0; j<NMU; j++) {
        QFtmp[j] = QF[j] - tmp2[j];
    }

    free(tmp1);
    free(tmp2);
}
