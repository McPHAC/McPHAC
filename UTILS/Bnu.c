#include "includes.h"


double Bnu(double localnu, double T) {
	// The blackbody function

	double  ans, arg;

	arg = HPOVERKB*localnu/T;
	ans = 2.0 * HPLANCK * pow(localnu,3.0)/SQR(CCC)/ (exp (arg) -1.0);
	return(ans);
}
