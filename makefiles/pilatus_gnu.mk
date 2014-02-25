CPP = mpicxx
GCC = mpicc

CFLAGS = -g -Wall -fopenmp #-dynamic
CXXFLAGS = -std=c++11 $(CFLAGS)

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /users/bruecks/GreenSolver/

# Common include paths
INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/QHULL/src/
INCAMD        = -I$(LIB_TOP)/AMD/Include/
INCMPS        = -I$(LIB_TOP)/MUMPS_4.10.0/include/
INCPDI        = -I$(LIB_TOP)/PDIV/
INCPOR        = -I$(LIB_TOP)/MUMPS_4.10.0/PORD/include/
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
METISLIB      = $(LIB_TOP)/Metis/libmetis.a
MUMPSLIB      = $(LIB_TOP)/MUMPS_4.10.0/lib/libzmumps.a
MUMPDLIB      = $(LIB_TOP)/MUMPS_4.10.0/lib/libdmumps.a
MUMPSCOM      = $(LIB_TOP)/MUMPS_4.10.0/lib/libmumps_common.a
PDIVLIB       = $(LIB_TOP)/PDIV/libpdiv.a
PORDLIB       = $(LIB_TOP)/MUMPS_4.10.0/lib/libpord.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a

# Directory of the rest of the libraries
LIB_TOP2      = /users/bruecks/

LINLIN = $(LIB_TOP2)/CSelInv/EXAMPLES/C2Finterface.o /users/bruecks/CSelInv/LIB/libcsupldlt.a

PARDISO_SO = $(LIB_TOP2)/Pardiso_SelInv/libpardiso491-GNU430-X86-64.so

DMALLOC = -L/apps/rosa/ddt/4.1.1/lib/64/ -ldmallocthcxx -z muldefs

LFLAGS = -L$(LIB_TOP2)/cp2k/lib/pilatus/popt/
DFLAGS = -DAdd_
LIBS = $(ARPACKLIB) -llapack -lblas $(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(MUMPDLIB) $(MUMPSCOM) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(QHULLLIB) \
	-lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lfftw3 -lm -lrt -fopenmp -lgfortran -lstdc++ 
INCLUDES = $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU)
