#SCOREP = scorep

CPP = $(SCOREP) CC
GCC = cc
NVCC = nvcc

#CFLAGS = -Wall -fopenmp -O3 -ffast-math -funroll-loops -march=native
#CFLAGS = -Wall -dynamic -fno-omit-frame-pointer -O1 -g -fsanitize=leak
CFLAGS = -Wall -g -fopenmp
CXXFLAGS = -std=c++11 $(CFLAGS)
NVCCFLAGS = -arch=compute_35 -code=sm_35

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s662/omendft_libraries

# Common include paths
INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/QHULL/src/
INCAMD        = -I$(LIB_TOP)/AMD/Include/
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include/
INCPEX        = -I$(LIB_TOP)/pexsi_v0.9.0/include/
INCMAG        = -I$(LIB_TOP)/magma-1.6.2p4_libsci_cuda7/include/

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a
PEXSILIB      = $(LIB_TOP)/pexsi_v0.9.0/src/libpexsi_daint.a
MAGMALIB      = $(LIB_TOP)/magma-1.6.2p4_libsci_cuda7/lib/libmagma.a

LINLIN = $(LIB_TOP)/CSelInv/EXAMPLES/C2Finterface.o $(LIB_TOP)/CSelInv/LIB/libcsupldlt.a

DMALLOC = -L/apps/common/ddt/6.0-Suse-11/lib/64/ -ldmallocthcxx -z muldefs

LFLAGS = -L$(LIB_TOP)/cp2k/lib/CRAY-XC30-gfortran-pexsi/popt/
DFLAGS = -DAdd_ -Dlibcp2k -DHAVE_MUMPS -DSPLITSOLVE
LIBS = $(UMFPACKLIB) $(AMDLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) $(PEXSILIB) $(MAGMALIB) \
	-lcp2k /project/ch5/alazzaro/libsmm/affinity/sandybridge_gcc_4.9.0/lib/libsmm_dnn_cray.gnu.a \
	-lzmumps -lsuperlu_dist -lptscotch -lptscotcherr -lscotch \
	-lcuda -lcudart -lcublas -lcufft -lcusparse \
	-lsci_gnu_mp -lm -lrt -fopenmp -lgfortran -lstdc++
INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCAZT) $(INCQHU) $(INCPEX) $(INCMAG)
