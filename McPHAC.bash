#!/bin/bash

# Physical Parameters
LOGTEFF=6.5   # Log Target effective temperature, K
GSURFACE=2.43e14   # Surface gravitational acceleration, cm s-2


# Computational parameters
#MINCOL=-9.0 # Min. log(y) considered for Teff<10^6.5K
#MAXCOL=2.0 # Initial max. log(y), increased if necessary to meet MAXTAUTAU criterion
MINCOL=-5.0 # Min. log(y) considered for Teff=10^6.5K
MAXCOL=2.0 # Initial max. log(y), increased if necessary to meet MAXTAUTAU criterion
#MINCOL=-8.0 # Min. log(y) considered in Z96 for Teff=10^5.3K
#MAXCOL=-0.6 # Max. log(y) considered in Z96 for Teff=10^5.3K (set MAXCOLTAU=0.01 to keep unchanged)
#MINCOL=-7.7 # Min. log(y) considered in Z96 for Teff=10^5.6K
#MAXCOL=0.2 # Max. log(y) considered in Z96 for Teff=10^5.6K (set MAXCOLTAU=0.01 to keep unchanged)
#MINCOL=-7.4 # Min. log(y) considered in Z96 for Teff=10^5.9K
#MAXCOL=0.9 # Max. log(y) considered in Z96 for Teff=10^5.9K (set MAXCOLTAU=0.01 to keep unchanged)
#MINCOL=-7.0 # Min. log(y) considered in Z96 for Teff=10^6.2K (not shown, so taken to be average of 5.9 and 6.5)
#MAXCOL=1.4 # Max. log(y) considered in Z96 for Teff=10^6.2K (set MAXCOLTAU=0.01 to keep unchanged)
#MINCOL=-6.6 # Min. log(y) considered in Z96 for Teff=10^6.5K
#MAXCOL=1.9 # Max. log(y) considered in Z96 for Teff=10^6.5K (set MAXCOLTAU=0.01 to keep unchanged)

MAXCOLTAU=80.0  # Consider log(y) large enough to have at least MAXCOLTAU optical depths at the largest freq.
#MAXCOLTAU=0.01  # Set this value of MAXCOLTAU to keep the largest y considered unchanged

TGUESSBDYCOND=0.264837817  # Ratio of T to Teff at the surface to use for initial temperature profile

NDEPTHS=200   # Initial number of depths points in the temperature correction
MAXFACTOR=1   # Maximum factor to multiply NDEPTHS by (NDEPTHS doubled until and including this factor)
NDEPTHSNU=200 # Initial number of depths points in the radiative transfer
MAXFACTORNU=1 # Maximum factor to multiply NDEPTHSNU by (NDEPTHSNU doubled until and including this factor)
NMU=10        # Number of mu points over range [0,1]
NFREQ=100     # Number of photon frequency bins

MAXFRACTEMPCHANGE=0.0001  # Continue iteration until max. fractional temp. change < MAXFRACTEMPCHANGE
MAXITER=20 # Maximum number of iterations allowed

ANIST=0  # Whether or not to treat Thomson scattering anisotropically (should only be set if FEAUTRIER is)

# Run McPHAC
cmnd="./McPHAC $LOGTEFF $GSURFACE $MINCOL $MAXCOL $NDEPTHS $MAXFACTOR $NDEPTHSNU $MAXFACTOR $NMU $NFREQ $MAXFRACTEMPCHANGE $MAXITER $ANIST $MAXCOLTAU $TGUESSBDYCOND"
echo $cmnd
$cmnd
