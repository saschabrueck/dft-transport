CPP = CC
GCC = cc

CFLAGS = -g -Wall -fopenmp
CXXFLAGS = -std=c++11 $(CFLAGS)

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /zhome/academic/HLRS/pri/iprsabru/
LIB_TOP2      = /zhome/academic/HLRS/pri/iprsabru/GreenSolver/

# Common include paths
INCAZT        = -I$(LIB_TOP2)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP2)/QHULL/src/
INCAMD        = -I$(LIB_TOP2)/AMD/Include/
INCUFC        = -I$(LIB_TOP2)/UFconfig/
INCUMF        = -I$(LIB_TOP2)/UMFPACK/Include
INCPEX        = -I$(LIB_TOP)/pexsi_v0.6.0/include

# Common library paths
ARPACKLIB     = $(LIB_TOP2)/ARPACK2/libarpack.a
AZTECLIB      = $(LIB_TOP2)/AZTEC/lib/libaztec.a
QHULLLIB      = $(LIB_TOP2)/QHULL/src/libqhull.a
AMDLIB        = $(LIB_TOP2)/AMD/Lib/libamd.a
UMFPACKLIB    = $(LIB_TOP2)/UMFPACK/Lib/libumfpack.a
PEXSI         = $(LIB_TOP)/pexsi_v0.6.0/src/libpexsi_hornet.a

LFLAGS = -L$(LIB_TOP)/cp2k/lib/CRAY-XC30-gfortran/psmp/
DFLAGS = -DAdd_
LIBS = $(UMFPACKLIB) $(AMDLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) $(PEXSI) \
	-lzmumps -lsuperlu_dist \
	-lcp2k \
	-lsci_gnu_mp -lm -lrt -fopenmp -lgfortran -lstdc++
INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCAZT) $(INCQHU) $(INCPEX)
