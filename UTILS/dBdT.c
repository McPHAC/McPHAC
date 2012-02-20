#include "includes.h"


double dBdT(double localnu,double T) {
	// Derivative of the blackbody function w.r.t. temperature

	double ans, arg, exparg, expargminusone, prefactor;

	prefactor=DBDTPREFACTOR;
	arg = HPOVERKB*localnu/T;

	exparg = exp(arg);
	expargminusone = (arg > 1e-4) ? exparg-1.0 : arg + 0.5*arg*arg;
	ans  = 2.0*prefactor* pow(localnu, 4)/(T*T) * ( (arg>100) ? 1./exparg : exparg/SQR(expargminusone) ); // exparg can overflow double
	return(ans);
}
