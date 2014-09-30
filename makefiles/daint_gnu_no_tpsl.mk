CPP = CC
GCC = cc

CFLAGS = -g -Wall -fopenmp -dynamic
CXXFLAGS = -std=c++11 $(CFLAGS) -DMKL_PARDISO  # or just any 'old' pardiso

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s503/omendft_libraries
LIB_TOP2      = /project/s503/OMEN_XC30

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

DMALLOC = -L/apps/rosa/ddt/4.1.1/lib/64/ -ldmallocthcxx -z muldefs

#PARDISO_SO = -L$(LIB_TOP)/Pardiso_SelInv -lpardiso491-GNU430-X86-64 -Wl,-rpath=$(LIB_TOP)/Pardiso_SelInv
PARDISO_SO = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LFLAGS = -L$(LIB_TOP)/cp2k1/lib/CRAY-XC30-gfortran/psmp/ -Wl,-rpath,/opt/fftw/3.3.4.0/sandybridge/lib/
DFLAGS = -DAdd_
LIBS = $(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(MUMPDLIB) $(MUMPSCOM) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(SUPERLULIB) $(QHULLLIB) $(ARPACKLIB) \
	-lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lsci_gnu_mp -lm -lrt -fopenmp -lgfortran -lstdc++ 
INCLUDES = $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU)
