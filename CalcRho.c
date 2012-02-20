#include "includes.h"

double CalcRho(double P, double T) {
    double ans, frac;
    float   t6, P12, R;
    void  opalrho_();  // Fortran program, so declare it here

    if (!OPALEOS || (USEMAXT && T>MAXT)) {
        ans = MPROTON* P/(KBOLTZMAN*T);   // The ideal gas law
        return(ans);
    } else {
	t6=T/1e6;  // Temperature should be in 1e6 K units for OPAL EOS
	P12=P/1e12; // Pressure should be in megabars, or 1e12 dynes/cm-2
	opalrho_(&P12, &t6,  &R);

	if (USEMAXT && T>0.9*MAXT) { // Smooth transition to ideal gas law
            frac = (T-0.9*MAXT)/MAXT;
	    ans= (1.0-frac)*R + frac*MPROTON* P/(KBOLTZMAN*T);
        } else {
	    ans=R;
        }
    }
    return(ans);
}
