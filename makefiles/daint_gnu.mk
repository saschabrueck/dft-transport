CPP = CC
GCC = cc

CFLAGS = -g -Wall -fopenmp
CXXFLAGS = -std=c++11 $(CFLAGS)

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s503/omendft_libraries
LIB_TOP2      = /project/s503/OMEN_XC30

# Common include paths
INCAZT        = -I$(LIB_TOP2)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP2)/QHULL/src/
INCAMD        = -I$(LIB_TOP2)/AMD/Include/
INCUFC        = -I$(LIB_TOP2)/UFconfig/
INCUMF        = -I$(LIB_TOP2)/UMFPACK/Include/

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AZTECLIB      = $(LIB_TOP2)/AZTEC/lib/libaztec.a
QHULLLIB      = $(LIB_TOP2)/QHULL/src/libqhull.a
AMDLIB        = $(LIB_TOP2)/AMD/Lib/libamd.a
UMFPACKLIB    = $(LIB_TOP2)/UMFPACK/Lib/libumfpack.a

LINLIN = $(LIB_TOP)/CSelInv/EXAMPLES/C2Finterface.o $(LIB_TOP)/CSelInv/LIB/libcsupldlt.a

DMALLOC = -L/apps/rosa/ddt/4.1.1/lib/64/ -ldmallocthcxx -z muldefs

LFLAGS = -L$(LIB_TOP)/cp2k1/lib/CRAY-XC30-gfortran/psmp/
DFLAGS = -DAdd_
LIBS = $(UMFPACKLIB) $(AMDLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) \
	-lzmumps -lsuperlu_dist \
	-lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lsci_gnu_mp -lm -lrt -fopenmp -lgfortran -lstdc++ 
INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCAZT) $(INCQHU)
