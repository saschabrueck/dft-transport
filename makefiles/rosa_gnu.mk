CPP = CC

CFLAGS = -g -Wall -DAdd_ 
CXXFLAGS = -std=c++11 $(CFLAGS)

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

LFLAGS = -L/users/bruecks/cp2k/lib/CRAY-XE6-gfortran-hwtopo/popt/ 
DFLAGS =
LIBS = $(SUPERLULIB) $(UMFPACKLIB) $(AMDLIB) $(METISLIB) $(MUMPSLIB) $(PORDLIB) \
	-lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lm -lgfortran -lstdc++	-lfftw3 -lsci_gnu 
INCLUDES = $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD)
