#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

// Physical Parameters (external variables set by input)
double  TEFF;     // Target effective temperature, K
double  GSURFACE; // Surface gravitational acceleration, cm s-2

// Computational parameters  (external variables set by input)
double  MINCOL;   // Min. log(y) considered for Teff=10^6.5K
double  MAXCOL;   // Initial max. log(y), increased if necessary to meet MAXTAUTAU criterion
int  NDEPTHS;     // Initial number of depths points in the temperature correction
int  MAXFACTOR;   // Maximum factor to multiply NDEPTHS by (NDEPTHS doubled until and including this factor)
int  NDEPTHSNU;   // Initial number of depths points in the radiative transfer
int  MAXFACTORNU; // Maximum factor to multiply NDEPTHSNU by (NDEPTHSNU doubled until and including this factor)
int  NMU;         // Number of mu points over range [0,1]
int  NFREQ;       // Number of photon frequency bins

double  MAXFRACTEMPCHANGE; // Continue iteration until max. fractional temp. change < MAXFRACTEMPCHANGE
int  MAXITER;           // Maximum number of iterations allowed


double  MAXCOLTAU;  // Consider log(y) large enough to have at least MAXCOLTAU optical depths at the largest freq.
#define  MAXTAU    80.0  // Consider MAXTAU optical depths at each frequency in the radiative transfer calculation
			 // (will not consider y larger than the largest value used in the temperature corr.)

#define  MINFREQKT  0.05   // Fraction of KBOLTZMAN/ERGSPERKEV*HZPER1KEVPHOTON*TEFF to use as minimum frequency
#define  MAXFREQKT  120.0  // Multiple of KBOLTZMAN/ERGSPERKEV*HZPER1KEVPHOTON*TEFF to use as maximum frequency

double  TGUESSBDYCOND;  // Ratio of T to Teff at the surface to use for initial temperature profile

#define  MINFREQINTERP  -2.0 // Minimum rescaled frequency to use when interpolating normalized spectra
#define  MAXFREQINTERP   2.0 // Maximum rescaled frequency to use when interpolating normalized spectra
#define  NFREQINTERP    1000 // Number of frequency points to interpolate normalized spectrum at


// Misc
#define  FILENAME 1000  // String length of a filename
#define SQR(X)   ((X)*(X))


// Filenames
#define  GFF_FILE  "./gffgu.dat"


// Options
int ANIST;  // Whether or not to treat Thomson scattering anisotropically (should only be set if FEAUTRIER is)


// Options that may be used in further development (some only partially implemented)
#define FREQKT  1   // Whether or not to give the max and min frequencies as fractions of k*Teff (fixed range otherwise)
#define USEBLAS  1  // Whether or not to use Basic Linear Algebra Subprograms for matrix multiplication
#define USELAPACK  1  // Whether or not to use Linear Algebra PACKage routines (needed for matrix inversion)
#define USEINTERPOL  1  // Whether or not to use different depth points in the radiative transfer calculation
			// for each frequency, and then interpolate to the temperature correction depth points
#define USELOGCOLNU  1  // Use logarithmic spacing in y for the radiative transfer calculation (if not set, linear spacing is used)
#define FEAUTRIER  1  // Whether or not to use the Feautrier solution for the radiative transfer (if not set, a Rybicki method is used)
#define AUER  1  // Whether or not to use the discretization scheme described in Auer 1976 (only implemented for Feautrier solution)
#define DIFFUSION  1  // Whether or not to use the diffusion approx. (currently P=B) as bottom bdy condition in Auer method
#define OPALEOS  1  // Whether or not to use the OPAL functions for the EOS (if not set, the ideal gas law is used)
#define THOMSON  1  // Whether or not to consider Thomson scattering (if not set, zero scattering opacity is used)
#define FREEFREE  1 // Whether or not to consider free-free scattering (if not set, zero free-free opacity is used)
#define GREY  0  // Whether or not to use a frequency independent free-free opacity in place of the regular opacity
#define DAMPDT  0  // Whether or not to damp the temperature correction by only applying a fraction of the computed correction
#define MINITERDAMP  2  // First iteration at which to consider DAMPDT
#define MAXTCHANGERATIO  0.5  // If the ratio of the present maxtempchange to the previous is larger than this value (and DAMPDT is set)
			      // only a fraction DTFRAC of the computed temperature correction will be applied in this and future iterations
#define DTFRAC  0.5  // Fraction of computed temperature correction to apply once MAXTCHANGERATIO criterion is met
#define MESSAGES  1  // Whether or not to print messages to stderr
#define OUTPUT  1  // Whether or not to output the radiative field to (big) files
#define USEMAXT    0  // Whether or not to use the ideal gas law for temperatures beyond MAXT


// Parameters pertaining to options that may be used in further development
#define  MINFREQ  15.4  // Minimum log(nu) to consider if FREQKT is not set (0.01 keV)
#define  MAXFREQ  18.4  // Maximum log(nu) to consider if FREQKT is not set (10 keV)
#define  MAXT    5e7  // If USEMAXT is set, the ideal gas law is used for temperatures larger than MAXT


// Constants
#define  ERGPEREV    1.60217646e-12 // erg/eV
#define  GGG        6.67259e-8   // Newton's Gravitational Constant, cgs
#define  CCC        2.997925e10  // speed of light, cm/s
#define  SIGMA      5.67051e-5    // Stefan-Boltzman constant, erg cm^-2 K^-4
#define  KBOLTZMAN  1.380658e-16  // erg/K
#define  HPLANCK    6.6260755e-27  // erg/s
#define  MPROTON    1.6726231e-24  // g
#define  MELECTRON  9.1093897e-28 // g
#define  SIGMA_T    6.65245e-25  // Thomson cross section, cm^2
#define  ECHARGE    4.80324e-10  // electron charge, cgs
#define  RY         13.6056923   // eV

#define  DBDTPREFACTOR 3.538222e-58 // HPLANCK^2/C^2/KBOLTZMAN
#define  HPOVERKB   4.799216e-11    // Ratio of the planck constant to the boltzman constant
#define  ERGSPERKEV  1.6e-9    // erg/keV
#define  HZPER1KEVPHOTON  2.4e17   // Frequency of a 1 keV photon

#define  LNTEN   2.302585093  // Natural logarithm of 10


// Struct definitions (some quantities may not be computed in present implementation)
typedef struct { // column_type struct for a specific depth
    // Quantities pertaining to the atmosphere structure
    double   logT; // Logarithm of temperature in K
    double   deltaT; // Computed temperature correction, K
    double   y; // Column depth y, g/cm^2
    double   logy; // Logarithm of y in g/cm^2
    double   dy;  // Change in y between present and next depth
    double   *tau; // Optical depth depth from surface, at each frequency
    double   *dtau;  // Change in optical depth between present and next depth, at each frequency
    double   rho; // Density, g/cm^3
    double   z; // Vertical position, cm
    double   pressure; // Pressure, ba

    // Quantities pertaining to the radiative transfer
    double   sigma;   // Scattering opacity
    double   kR;  // Rosseland mean opacity
    double   *kappa;  // Free-free opacity, for each frequency
    double   *k; // Total opacity, as a function of frequency
    double   *rho_opacity; // Dimensionless scattering albedo
    double   **I;  // Specific intensity, at each angle and frequency
    double   **JanisT;  // Weighted mean intensity for anisotropic Thomson scattering, for each angle and freq.
    double   **u;  // u = P = 0.5 (I+ + I-), for each angle and frequency
    double   *S;  // Source function
    double   *J;  // Mean intensity
    double   *Jt;  // Mean intensity from the temperature correction
    double   *f; // Eddington factor f
    double   *h; // Eddington factor h
    double   *q; // Angle-averaged ratio of JanitT to regular J (should = 1 everywhere)
    double   *F;  // Flux, at each frequency
    double   totalFlux;   // Total flux
    double   Teff; // Effective temperature

    double   Jbar; // Quantity used in the Rybicki solution to the radiative transfer
    double   Jtbar; // Quantity used in the Rybicki solution to the temperature correction
    double   Jbarb; // Quantity used in the Rybicki solution to the radiative transfer, possibly the same as Jbar

    double   **BFtmp; // Matrix used in the Feautrier solution to the radiative transfer
    double   **BFtmpinv; // Matrix used in the Feautrier solution to the radiative transfer
    double   *QFtmp; // Vector used in the Feautrier solution to the radiative transfer
} column_type ;

typedef struct { // matrix_type struct containing matrices used for calculations
    double   **T; // Matrix used for Rybicki type solutions
    double   **TInverse; // Matrix used for Rybicki type solutions
    double   **U; // Matrix used for Rybicki type solutions
    double   **V; // Matrix used for Rybicki type solutions
    double    *K; // Vector used for Rybicki type solutions
    double   **W; // Matrix used for Rybicki type solutions
    double   *Q;  // Vector used for Rybicki type solutions

    double   **AF; // Matrix used for Feautrier type solutions
    double   **BF; // Matrix used for Feautrier type solutions
    double   **CF; // Matrix used for Feautrier type solutions
    double   *QF;  // Vector used for Feautrier type solutions
} matrix_type;


// Function declarations
column_type * AllocateColumns(int, int, int);
void     AllocateMatrix(matrix_type*, int, int);
double   Kappa_nu(double, double, double);
double   Bnu(double, double);

void    CalcJbar(double **, double *, column_type  *, int);
void    CalcK(double *, int, int, int, column_type *, int);
void    CalcU(double **, int, column_type *, int);
void    CalcT(double **, int, int, column_type *, int);
void    CalcV(double **, int, int, column_type *, int);

void    CalcAF(double **, int, int, column_type *, int);
void    CalcBF(double **, int, int, column_type *, int);
void    CalcBFtmp(double **, double **, double **, double **, double **);
void    CalcCF(double **, int, int, column_type *, int);
void    CalcQF(double *, int, int, int, column_type *, int);
void    CalcQFtmp(double *, double **, double **, double *, double *);
void    UpdateuF(int, matrix_type, column_type *, int);

void    CalcJtbar(double **, double *, column_type  *, int);
void    CalcKt(double *, int, int, column_type *, int);
void    CalcUt(double **, int, column_type *, int);
void    CalcTt(double **, int, int, column_type *, int);
void    CalcVt(double **, int, int, column_type *, int);
void    CalcEddingtonFactors(column_type*, double*, double*, int, int, int);
void    InterPolEdd(column_type *, int, int, column_type*, int, int);
void    OutputEddingtonFactors(column_type*, int, int, double, int, int);
void    UpdateJt(int, double**, double*, double**, column_type*, int, int, int);

void     CalcFluxes(column_type* , double*, double*, int, int, int);
void     CalcFluxest(column_type* , double*, double*, int, int);
void     CalcTeff(column_type  *, double *, double *, int, int);
void     CalcFluxesOld(column_type  *, int, int);
void     CalckR(column_type  *);
void     CalcMeanIntensityJ(column_type  *, double *, int, int, int);
void     CalcOpacities(column_type  *, int, int, int);
double   GauntFactor(double , double);
double   CalcRho(double, double);
void     CalcSpecificIntensities(column_type *, double *, int, int, int);
void     CalcSourceFunction(column_type *);
double   DeltaT(column_type *, int, int, double, int *);
void     PrintInterpolT(column_type *, int, column_type *, int, int);
double   dBdT(double,double);
double   * dvector(int);
double   ** ddvector(int);
void     free_vector(double*, long, long);
void    FreeColumns(column_type*, int, int, int);
void     FreeMatrix(matrix_type*, int, int);
column_type* GetColumnsLog(column_type*, int, int, int, int*, int*, int*, int);
column_type* GetColumnsNu(column_type *, int, int, int, int, int *);
void     InitW(double **, int);
void     InterPolJ(column_type *, int, double *, int, int, column_type *, int);
void     InterPolF(column_type *, int, double *, int, int, column_type *, int);
void     InvertMatrix(double**, double**, int);
void     InvertTridiagonalMatrix(double**, double**, int);
void     MatrixMultiply(double**, double**, double**, int);
void     MatrixDiagonalMatrixMultiply(double**, double**, double**, int);
void     DiagonalMatrixMatrixMultiply(double**, double**, double**, int);
void     MatrixVectorMultiply(double *, double **, double *, int);
void     OutputSpectrum(column_type*, FILE*, int, double*, int, double*, int);
void     OutputSpectrumOld(column_type*, int);
void     OutputFluxes(column_type  *, double *, double *, int, int, int, int);
void     interpolate(double*, double*, int, double, double*, double*);
void     InterpolateArray(double*, double*, int, double*, double*, int, int, double, int, double);
void     cspline_set(int32_t, double*, double*, int32_t, double, int32_t, double, double*);
void     cspline_val(int32_t, double*, double*, double*, double, double*, double*, double*);
void     PrintOutTau(column_type*);
void     PrintOutOpacities(column_type*);
double   Sigma(double);
void      UpdateW(double**, double**, double**, double**, int);
void      UpdateQ(double *, double **, double **, double *, int);
void      Updateu(int, int, double**, double*, double**, column_type*, int);


// External variables (should perhaps be phased out)
double  *nu, *lognu, *dnu;
double  *mu, *dmu;
double  *dmudnu;
int     NDepths;
matrix_type  M;
column_type   *cfori;
column_type   cfornu;
double **gauntff, **gauntfactor, *ugrid, *g2grid;
int NGaunt, NGauntu, NGauntg2;
double  extern_maxtemperaturechange;
