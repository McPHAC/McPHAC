####################
# This is the Makefile for McPHAC
####################


####################
# Define directories where various files are stored

# BIN is the directory where the executable "McPHAC" will be stored
BIN  = .

# NS is the directory containing the McPHAC source code
NS   = .

# OPAL is the directory containing the OPAL Equation of State code
OPAL = ./OPAL

# UTILS is the directory where simple utility programs are stored
UTILS = ./UTILS

# LAPACK  contains the source code for select LAPACK routines
LAPACK = ./LAPACK

# BLAS contains the source code for the BLAS routines
BLAS = ./BLAS

# OPTIMIZEDLIBS contains pre-compiled optimized numerical libraries
OPTIMIZEDLIBS = ./OPTIMIZEDLIBS


####################
# Set various options

# Enable warnings
#WARNINGS = -Wall

# Enable gdb debugging info
DEBUG    = -ggdb

# To compile numerical routines from source, use the following LIB definition:
LIB  	 = -lgfortran -lm -L $(BLAS)/ -lblas
# To use pre-compilied optimized AMD numerical routines, use the following LIB definition:
#LIB  	  = -lgfortran -lm -L $(OPTIMIZEDLIBS)/ -lacml

# Store profiling info
#PROFILE  = -pg

# Enable optimization
OPTIMIZE = -O3

# Set compiler flags
CFLAGS   = $(WARNINGS) $(DEBUG) $(PROFILE) $(OPTIMIZE) -I.
FFLAGS   = $(WARNINGS) $(DEBUG) $(PROFILE) $(OPTIMIZE)
F90FLAGS = $(WARNINGS) $(DEBUG) $(PROFILE) $(OPTIMIZE)

# Set compilers
CC      =      gcc
FORTRAN =      gfortran
FORTRAN90 =    gfortran

# Set linker flags
LDFLAGS =      $(LIB) $(DEBUG) $(PROFILE)


####################
# Group related object files

# General utilities
UTILOBJS = $(UTILS)/dvector.o \
	$(UTILS)/ddvector.o \
	$(UTILS)/MatrixMultiply.o \
	$(UTILS)/MatrixDiagonalMatrixMultiply.o \
	$(UTILS)/DiagonalMatrixMatrixMultiply.o \
	$(UTILS)/MatrixVectorMultiply.o \
	$(UTILS)/InvertMatrix.o \
	$(UTILS)/InvertTridiagonalMatrix.o \
	$(UTILS)/interpolate.o \
	$(UTILS)/InterpolateArray.o \
	$(UTILS)/cspline_set.o \
	$(UTILS)/cspline_val.o \
	$(UTILS)/spline.o \
	$(UTILS)/quadrule.o \
	$(UTILS)/AllocateMatrix.o \
	$(UTILS)/FreeMatrix.o \
	$(UTILS)/FreeColumns.o \
	$(UTILS)/Bnu.o \
	$(UTILS)/dBdT.o \
	$(UTILS)/AllocateColumns.o  

# Lapack routines
LAPACKOBJ =   $(LAPACK)/dgesv.o \
	$(LAPACK)/dgetf2.o \
	$(LAPACK)/dgetrf.o \
	$(LAPACK)/dgetrs.o \
	$(LAPACK)/dgtsv.o \
	$(LAPACK)/dlamch.o \
	$(LAPACK)/ilaenv.o \
	$(LAPACK)/ieeeck.o \
	$(LAPACK)/iparmq.o \
	$(LAPACK)/dgetri.o \
	$(LAPACK)/dtrti2.o \
	$(LAPACK)/dtrtri.o \
	$(LAPACK)/dlaswp.o

# McPHAC functions
NSOBJS 	= $(OPAL)/OPALrho.o \
	$(OPAL)/IEOSxtrin_H_He.o \
	$(NS)/GetColumnsLog.o \
	$(NS)/GetColumnsNu.o \
	$(NS)/GauntFactor.o \
	$(NS)/CalcAF.o \
	$(NS)/CalcBF.o \
	$(NS)/CalcBFtmp.o \
	$(NS)/CalcCF.o \
	$(NS)/CalcQF.o \
	$(NS)/CalcQFtmp.o \
	$(NS)/UpdateuF.o \
	$(NS)/InitW.o \
	$(NS)/CalcU.o \
	$(NS)/CalcUt.o \
	$(NS)/CalcT.o \
	$(NS)/CalcTt.o \
	$(NS)/CalcV.o \
	$(NS)/CalcVt.o \
	$(NS)/CalcK.o \
	$(NS)/CalcKt.o \
	$(NS)/UpdateW.o \
	$(NS)/UpdateQ.o \
	$(NS)/CalcJbar.o \
	$(NS)/CalcJtbar.o \
	$(NS)/Updateu.o \
	$(NS)/CalcMeanIntensityJ.o \
	$(NS)/CalcEddingtonFactors.o \
	$(NS)/InterPolEdd.o \
	$(NS)/OutputEddingtonFactors.o \
	$(NS)/InterPolJ.o \
	$(NS)/CalcFluxes.o \
	$(NS)/OutputFluxes.o \
	$(NS)/CalcFluxest.o \
	$(NS)/InterPolF.o \
	$(NS)/CalcTeff.o \
	$(NS)/CalcFluxesOld.o \
	$(NS)/OutputSpectrum.o \
	$(NS)/OutputSpectrumOld.o \
	$(NS)/UpdateJt.o \
	$(NS)/DeltaT.o \
	$(NS)/Kappa_nu.o \
	$(NS)/Sigma.o \
	$(NS)/CalcOpacities.o \
	$(NS)/CalckR.o \
	$(NS)/CalcRho.o \
	$(NS)/PrintInterpolT.o \
	$(NS)/PrintOutTau.o \
	$(NS)/PrintOutOpacities.o \
	$(NS)/McPHAC.o


# Matrix inversion testing routines
INVOBJS = $(UTILOBJS) $(UTILS)/invert.o


####################
# Define make rules

# If no arguments or argument 'all', make McPHAC and invert
all:    McPHAC invert

# Use this rule for McPHAC if numerical routines are compiled from source 
McPHAC	: $(NSOBJS)  $(LAPACKOBJ) $(UTILOBJS) blaslib
	$(CC) $(NSOBJS) $(UTILOBJS) $(LAPACKOBJ) $(LDFLAGS) -o $(BIN)/$@
# Use this rule for McPHAC if pre-compiled optimized numerical routines are used
#McPHAC	: $(NSOBJS) $(UTILOBJS)
#	$(CC) $(NSOBJS) $(UTILOBJS) $(LDFLAGS) -o $(BIN)/$@

# Define rule to compile BLAS routines
blaslib: 
	( cd $(BLAS); $(MAKE) double)

# Set dependency on includes.h, where all parameters are defined
$(NSOBJS) $(UTILOBJS) $(INVOBJS) : $(NS)/includes.h

# Define implicit rules for making object files
.f.o:
	$(FORTRAN) $(FFLAGS) -c $< -o $@

.SUFFIXES: .o .f90
.f90.o:
	$(FORTRAN90) $(F90FLAGS) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

# Use this rule for the matrix inversion tester if compiling numerical routines from source
invert	: $(INVOBJS) $(LAPACKOBJ) blaslib
	 $(CC) $(LAPACKOBJ) $(INVOBJS) $(LIB) $(LDFLAGS) -o $(BIN)/$@
# Use this rule for the matrix inversion tester if using pre-compiled numerical routines
#invert	: $(INVOBJS)
#	$(CC) $(INVOBJS) $(LIB) $(LDFLAGS) -o $(BIN)/$@

# Define rule for cleaning directory (leading dash prevents make from exiting on rm error)
clean	: 
	-rm McPHAC invert
	-rm $(BIN)/*.o $(BIN)/*~
	-rm $(NS)/*.o $(NS)/*~
	-rm $(OPAL)/*.o $(OPAL)/*~
	-rm $(UTILS)/*.o $(UTILS)/*~
	-rm $(LAPACK)/*.o $(LAPACK)/*~
	-rm $(BLAS)/*.o $(BLAS)/*~ $(BLAS)/*.a
