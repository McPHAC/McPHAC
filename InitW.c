#include "includes.h"

void InitW(double **W, int n) {
    // Initializes the W matrix to the negative identity matrix

    int i;

    for (i=0; i<n; i++) {
        memset(W[i], 0, n*sizeof(double));
        W[i][i]=-1.0;
    }
}

