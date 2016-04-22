CPP = CC
GCC = cc

#CFLAGS = -Wall -fopenmp -O3 -ffast-math -funroll-loops -march=native
CFLAGS = -Wall -fopenmp -g
CXXFLAGS = -std=c++11 $(CFLAGS) -DMKL_PARDISO  # or just any 'old' pardiso

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s579/omendft_libraries
LIB_TOP2      = /project/s579/OMEN_XC30

# Common include paths
INCAZT        = -I$(LIB_TOP2)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP2)/QHULL/src/
INCAMD        = -I$(LIB_TOP2)/AMD/Include/
INCSLU        = -I$(LIB_TOP2)/SuperLU_DIST_2.0/SRC/
INCMPS        = -I$(LIB_TOP2)/MUMPS_4.10.0/include/
INCPDI        = -I$(LIB_TOP2)/PDIV/
INCPOR        = -I$(LIB_TOP2)/MUMPS_4.10.0/PORD/include/
INCSLU        = -I$(LIB_TOP2)/SuperLU_DIST_2.0/SRC
INCUFC        = -I$(LIB_TOP2)/UFconfig/
INCUMF        = -I$(LIB_TOP2)/UMFPACK/Include

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AMDLIB        = $(LIB_TOP2)/AMD/Lib/libamd.a
AZTECLIB      = $(LIB_TOP2)/AZTEC/lib/libaztec.a
METISLIB      = $(LIB_TOP2)/Metis/libmetis.a
MUMPSLIB      = $(LIB_TOP2)/MUMPS_4.10.0/lib/libzmumps.a
MUMPDLIB      = $(LIB_TOP2)/MUMPS_4.10.0/lib/libdmumps.a
MUMPSCOM      = $(LIB_TOP2)/MUMPS_4.10.0/lib/libmumps_common.a
PDIVLIB       = $(LIB_TOP2)/PDIV/libpdiv.a
PORDLIB       = $(LIB_TOP2)/MUMPS_4.10.0/lib/libpord.a
QHULLLIB      = $(LIB_TOP2)/QHULL/src/libqhull.a
UMFPACKLIB    = $(LIB_TOP2)/UMFPACK/Lib/libumfpack.a
SUPERLULIB    = $(LIB_TOP2)/SuperLU_DIST_2.0/Lib/libsuperlu_dist_2.0.a

LINLIN = $(LIB_TOP)/CSelInv/EXAMPLES/C2Finterface.o $(LIB_TOP)/CSelInv/LIB/libcsupldlt.a

PARDISO_SO = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LFLAGS = -L$(LIB_TOP)/cp2k/lib/CRAY-XC30-gfortran/popt/
DFLAGS = -DAdd_ -DHAVE_SUPERLU -DHAVE_UMFPACK -Dlibcp2k
LIBS = $(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(MUMPSCOM) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(SUPERLULIB) $(QHULLLIB) $(ARPACKLIB) \
	-lcp2k /project/ch5/vondele/libsmm_alfio/affinity/sandybridge_gcc_4.9.0/lib/libsmm_dnn_cray.gnu.a \
	-lsci_gnu_mp -lm -lrt -fopenmp -lgfortran -lstdc++ 
INCLUDES = $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU)
