# module switch PrgEnv-pgi/5.2.82 PrgEnv-gnu
# module load cray-tpsl/16.03.1
# module load cudatoolkit
# module load fftw

#SCOREP = scorep

CPP = $(SCOREP) CC
GCC = cc
NVCC = nvcc

CFLAGS = -Wall -O3 -ffast-math -funroll-loops -march=native
#CFLAGS = -Wall -dynamic -fno-omit-frame-pointer -O1 -g -fsanitize=leak
#CFLAGS = -Wall -g -fopenmp
CXXFLAGS = -std=c++11 $(CFLAGS)
NVCCFLAGS = -arch=compute_35 -code=sm_35

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_SAB       = /ccs/home/luisier/libraries
LIB_TOP       = /ccs/home/luisier/OMEN_XK7

# Common include paths
INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/QHULL/src/
INCAMD        = -I$(LIB_TOP)/AMD/Include/
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include/
INCMAG        = -I$(LIB_TOP)/magma-1.6.2p5/include
INCPEX        = -I$(LIB_SAB)/pexsi_v0.9.2/include/

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK/libarpack.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a
MAGMALIB      = $(LIB_TOP)/magma-1.6.2p5/lib/libmagma.a
PEXSILIB      = $(LIB_SAB)/pexsi_v0.9.2/src/libpexsi_titan.a

LFLAGS = -L$(LIB_SAB)/cp2k/lib/CRAY-XK7-gfortran-pexsi/popt/
DFLAGS = -DAdd_ -Dlibcp2k -DHAVE_MUMPS -DHAVE_SUPERLU -DHAVE_PEXSI -DHAVE_SPLITSOLVE -DHAVE_OMEN_POISSON
LIBS = -lcp2k \
	$(UMFPACKLIB) $(AMDLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) $(PEXSILIB) $(MAGMALIB) \
	-lzmumps -lsuperlu_dist -lptscotch -lptscotcherr -lscotch \
	-lcuda -lcudart -lcublas -lcufft -lcusparse \
	-lsci_gnu -lm -lrt -fopenmp -lgfortran -lstdc++
INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCAZT) $(INCQHU) $(INCPEX) $(INCMAG)
