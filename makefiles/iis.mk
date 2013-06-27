CPP = mpicxx

CFLAGS = -g -Wall -DAdd_ 
CXXFLAGS = -std=c++11 $(CFLAGS)

BASE     = /home/mauro/sw
GCC_LIB  = $(BASE)/gcc/4.8.1/lib64
SCALAPACK_LIB = $(BASE)/scalapack/lib
OPENBLAS_LIB = $(BASE)/openblas/lib
FFTW_LIB = $(BASE)/fftw/gcc/lib
FFTW_INC = $(BASE)/fftw/gcc/include

# Directory of the libraries used by OMEN
LIB_TOP       = /usr/zupo/sim3/bruecks/GreenSolver/

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

LFLAGS = -L/usr/zupo/sim3/bruecks/cp2k/cp2k/lib/iis_ee-x86-64-gcc-openblas-scalapack/popt/ \
	-L$(GCC_LIB) -L$(SCALAPACK_LIB) -L$(OPENBLAS_LIB) -L$(FFTW_LIB)
DFLAGS =
LIBS = $(SUPERLULIB) $(UMFPACKLIB) $(AMDLIB) $(METISLIB) $(MUMPSLIB) $(PORDLIB) \
	-lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lm -lgfortran -lstdc++	-lfftw3 -lscalapack -lopenblas 
INCLUDES = $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) -I$(FFTW_INC)
