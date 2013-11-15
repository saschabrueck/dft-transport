CPP = CC
GCC = cc

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
INCSLU        = -I$(LIB_TOP)/SuperLU_DIST_2.0/SRC/
INCMPS        = -I$(LIB_TOP)/MUMPS/include/
INCPDI        = -I$(LIB_TOP)/PDIV/
INCPOR        = -I$(LIB_TOP)/MUMPS/PORD/include/
INCSLU        = -I$(LIB_TOP)/SuperLU_DIST_2.0/SRC
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK/libarpack.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
METISLIB      = $(LIB_TOP)/Metis/libmetis.a
MUMPSLIB      = $(LIB_TOP)/MUMPS/lib/libzmumps.a
PDIVLIB       = $(LIB_TOP)/PDIV/libpdiv.a
PORDLIB       = $(LIB_TOP)/MUMPS/lib/libpord.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
SUPERLULIB    = $(LIB_TOP)/SuperLU_DIST_2.0/Lib/libsuperlu_dist_2.0.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a

LINLIN = /users/bruecks/CSelInv/EXAMPLES/C2Finterface.o /users/bruecks/CSelInv/LIB/libcsupldlt.a

PARDISO_SO = /users/bruecks/bin/libpardiso491-GNU430-X86-64.so

LFLAGS = -L/users/bruecks/cp2k/lib/CRAY-XE6-gfortran-hwtopo/popt/ #-Wl,-rpath,/opt/cray/mpt/5.6.1/gni/mpich2-gnu/47/lib/ -Wl,-rpath,/opt/fftw/3.3.0.1/interlagos/lib/
#LFLAGS = -L/users/bruecks/cp2k2/cp2k/cp2k/lib/CRAY-XE6-gfortran-hwtopo/popt/
DFLAGS = -DAdd_
LIBS = $(SUPERLULIB) $(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(QHULLLIB) \
	-lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lsci_gnu_mp -lm -lrt -fopenmp -lgfortran -lstdc++ 
INCLUDES = $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU)
