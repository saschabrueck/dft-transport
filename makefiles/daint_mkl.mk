# module load daint-gpu
# module switch PrgEnv-cray PrgEnv-gnu
# module load fftw
# module load ddt
# module load cray-tpsl
# module load cudatoolkit

#SCOREP = scorep

CPP = $(SCOREP) CC
GCC = cc
NVCC = nvcc

CFLAGS = -Wall -fopenmp -O3 -ffast-math -funroll-loops -march=native
#CFLAGS = -Wall -dynamic -fno-omit-frame-pointer -O1 -g -fsanitize=leak
#CFLAGS = -Wall -g -fopenmp
CXXFLAGS = -std=c++11 $(CFLAGS)
NVCCFLAGS = -arch=compute_60 -code=sm_60

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s662/omendft_libraries

# Common include paths
INCAZT        = -I$(LIB_TOP)/../OMEN_XC50/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/../OMEN_XC50/QHULL/src/
INCAMD        = -I$(LIB_TOP)/../OMEN_XC50/AMD/Include/
INCUFC        = -I$(LIB_TOP)/../OMEN_XC50/UFconfig/
INCUMF        = -I$(LIB_TOP)/../OMEN_XC50/UMFPACK/Include/
INCPEX        = -I$(LIB_TOP)/pexsi_v0.10.1/include/

# Common library paths
AZTECLIB      = $(LIB_TOP)/../OMEN_XC50/AZTEC/lib/libaztec.a
QHULLLIB      = $(LIB_TOP)/../OMEN_XC50/QHULL/src/libqhull.a
AMDLIB        = $(LIB_TOP)/../OMEN_XC50/AMD/Lib/libamd.a
UMFPACKLIB    = $(LIB_TOP)/../OMEN_XC50/UMFPACK/Lib/libumfpack.a
PEXSILIB      = $(LIB_TOP)/pexsi_v0.10.1/src/libpexsi_daint.a

DMALLOC = -L/apps/common/UES/SLES12/ddt/6.1.2/lib/64/ -ldmallocthcxx -z muldefs

LFLAGS = -L$(LIB_TOP)/cp2k/lib/pexsi/popt/
DFLAGS = -DAdd_ -Dlibcp2k -DHAVE_MUMPS -DHAVE_SUPERLU -DHAVE_PEXSI -DHAVE_SPLITSOLVE -DHAVE_OMEN_POISSON
LIBS = -lcp2k \
	$(UMFPACKLIB) $(AMDLIB) $(AZTECLIB) $(QHULLLIB) $(PEXSILIB) \
	-lzmumps -lsuperlu_dist -lptscotch -lptscotcherr -lscotch -lmagma\
	-lcuda -lcudart -lcublas -lcufft -lcusparse \
	-lmkl_gf_lp64 -lmkl_sequential -lmkl_core \
	-lm -lrt -fopenmp -lgfortran -lstdc++
INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCAZT) $(INCQHU) $(INCPEX) $(INCMAG)
