#include "includes.h"
#define N 500

#define ERRORSONLY 1
#define TRIDIAG 0
#define FRACERR 1e-9
int main(void) {

    double   **in, **out, **product;
    int   i, j;

    in = (double **) malloc (N*sizeof(double * ));
    out = (double **) malloc (N*sizeof(double * ));
    product = (double **) malloc (N*sizeof(double * ));
    for (i=0; i<N; i++) {
        in[i] = (double *)malloc(N*sizeof(double));
        out[i] = (double *)malloc(N*sizeof(double));
        product[i] = (double *)malloc(N*sizeof(double));
        for (j=0; j<N; j++) {
            if (TRIDIAG) {
		in[i][j] = (abs(i-j)<2) ? (((float)(rand())/(float)(RAND_MAX))-0.5) : 0.0;
            } else {
	        in[i][j] = (((float)(rand())/(float)(RAND_MAX))-0.5);
            }
            if (!(ERRORSONLY)) fprintf(stdout,"%f\t", in[i][j]);
        }
        if (!(ERRORSONLY)) fprintf(stdout,"\n");
    }


    if (TRIDIAG) {
        InvertTridiagonalMatrix(out, in, N);
    } else {
        InvertMatrix(out, in, N);
    }

    if (!ERRORSONLY) {
        fprintf(stdout,"\n\n");

        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                fprintf(stdout,"%f\t", out[i][j]);
            }
            fprintf(stdout,"\n");
        }
    }

    MatrixMultiply(product,in,out,N);

    fprintf(stdout,"\n\n");
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (!(ERRORSONLY))  fprintf(stdout,"%f\t", product[i][j]);
            if ((i==j && fabs(product[i][j]-1.0) > FRACERR) || (i!=j && fabs(product[i][j])>FRACERR)) {
                fprintf(stdout,"***** ERRROR ****");
            }
        }
        if (!(ERRORSONLY))  fprintf(stdout,"\n");
    }

return 0;
}
