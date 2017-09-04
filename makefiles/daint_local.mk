# module load daint-gpu; module unload PrgEnv-cray; module load PrgEnv-gnu/6.0.3; module load fftw/3.3.4.10; module load cudatoolkit; module load magma; module load ddt; module load cray-tpsl

LEX  = flex
YACC = bison

#SCOREP = scorep
CPP  = $(SCOREP) CC
GCC  = cc
NVCC = nvcc

CFLAGS    = -Wall -fopenmp -O3 -ffast-math -funroll-loops -march=native
CXXFLAGS  = -std=c++11 $(CFLAGS)
NVCCFLAGS = -arch=compute_60 -code=sm_60

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s662/omendft_libraries
LDIR          = /project/s662/hbani/cp2k-omen/cp2k/cp2k/tools/toolchain/install

# Common include paths
INCHYP        = -I/project/s662/hbani/cp2k-omen/omen/install/libs/hypre/include/
INCQHU        = -I$(LIB_TOP)/../OMEN_XC50/QHULL/src/
INCAMD        = -I$(LIB_TOP)/../OMEN_XC50/AMD/Include/
INCUFC        = -I$(LIB_TOP)/../OMEN_XC50/UFconfig/
INCUMF        = -I$(LIB_TOP)/../OMEN_XC50/UMFPACK/Include/
INCMAG        = -I$(LIB_TOP)/magma-2.2.0/include/
INCPEX        = -I$(LDIR)/pexsi-0.10.1/include/

# Common library paths
HYPRELIB      = -L/project/s662/hbani/cp2k-omen/omen/install/libs/hypre/lib/
QHULLLIB      = -L$(LIB_TOP)/../OMEN_XC50/QHULL/src/ 
UMFPACKLIB    = -L$(LIB_TOP)/../OMEN_XC50/UMFPACK/Lib/ 
AMDLIB        = -L$(LIB_TOP)/../OMEN_XC50/AMD/Lib/ 
MAGMALIB      = -L$(LIB_TOP)/magma-2.2.0/lib/ 
PEXSILIB      = -L$(LDIR)/pexsi-0.10.1/lib/
LIBELPA       = -L$(LDIR)/elpa-2016.05.004/lib/
LIBXSMM       = -L$(LDIR)/libxsmm-1.6.4/lib/
LIBXC         = -L$(LDIR)/libxc-2.2.2/lib/ 
LIBINT        = -L$(LDIR)/libint-1.1.4/lib/ 

# CP2K library path
CP2KLIB  = -L/project/s662/hbani/cp2k-omen/cp2k/cp2k/lib/local_new/psmp

DMALLOC  = -L/apps/common/UES/SLES12/ddt/6.1.2/lib/64/ -ldmallocthcxx -z muldefs

LFLAGS   = $(CP2KLIB) $(LIBELPA) $(LIBXSMM) $(LIBINT) $(LIBXC) $(UMFPACKLIB) $(AMDLIB) $(HYPRELIB) $(QHULLLIB) $(PEXSILIB) $(MAGMALIB)
DFLAGS   = -DAdd_  

LIBS  = -lcp2k -lelpa_openmp -lxsmmf -lxsmm -ldl -lxcf90 -lxc -lderiv -lint -lpexsi \
	-lHYPRE -lqhull -lumfpack -lamd \
	-lzmumps -lsuperlu_dist -lptscotch -lptscotcherr -lscotch \
	-lcuda -lcudart -lcublas -lcufft -lcusparse -lmagma \
	-lm -lrt -fopenmp -lgfortran -lstdc++

INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCHYP) $(INCQHU) $(INCPEX) $(INCMAG)

