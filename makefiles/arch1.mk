# source  /home/seyedb/cp2k-omen/cp2k/cp2k/tools/toolchain/install/setup
#
LEX      = flex
YACC     = bison

CXX      = mpic++ 
CC       = gcc 
NVCC     = /usr/local/cuda-7.5/bin/nvcc 

CFLAGS    = -g -w -Wall -fopenmp 
CXXFLAGS  = -std=c++11 $(CFLAGS)
NVCCFLAGS = -Xcompiler=--std=gnu++98 -D__GNUC__=4 -D__GNUC_MINOR__=9 -w

TOP_DIR       = /home/seyedb/cp2k-omen
LIB_TOP       = $(TOP_DIR)/libs
TOOLCHAIN     = /home/seyedb/cp2k-omen/cp2k/cp2k/tools/toolchain

# Common include paths
INCSSPARSE    = $(LIB_TOP)/SuiteSparse/include/
INCMUMPS      = $(LIB_TOP)/MUMPS/include/
INCHYPRE      = $(LIB_TOP)/hypre/include/
INCQHULL      = $(LIB_TOP)/qhull/include/libqhull/
INCMAGMA      = $(LIB_TOP)/magma/include/
INCSLUDIST    = $(TOOLCHAIN)/install/superlu_dist-5.1.2/include/
INCPEXSI      = $(TOOLCHAIN)/install/pexsi-0.10.1/include/

INCCUDA       = /usr/local/cuda-7.5/include/

# Common library paths
LIBSSPARSE    = $(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig 
LIBMUMPS      = $(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack 
LIBHYPRE      = $(LIB_TOP)/hypre/lib -lHYPRE 
LIBQHULL      = $(LIB_TOP)/qhull/lib -lqhullstatic 
LIBMAGMA      = $(LIB_TOP)/magma/lib -lmagma
LIBPARMETIS   = $(TOOLCHAIN)/install/parmetis-4.0.3/lib -lparmetis -lmetis 
LIBSLUDIST    = $(TOOLCHAIN)/install/superlu_dist-5.1.2/lib -lsuperlu_dist 
LIBPEXSI      = $(TOOLCHAIN)/install/pexsi-0.10.1/lib -lpexsi 
LIBCP2K       = /home/seyedb/cp2k-omen/cp2k/cp2k/lib/local/popt -lcp2k -lxsmmf -lxsmm -lderiv -lint -lxcf90 -lxc -lfftw3 

LIBCUDA       = /usr/local/cuda-7.5/lib64/ -lcudart -lcublas -lcusparse -lblas

LFLAGS   = -L$(LIBCP2K) -L$(LIBPEXSI) -L$(LIBPARMETIS) -L$(LIBSLUDIST) -L$(LIBSSPARSE) -L$(LIBMUMPS) -L$(LIBHYPRE) -L$(LIBQHULL) -L$(LIBCUDA) -L$(LIBMAGMA)

DFLAGS   = -DAdd_ 

LIBS     = -lrt -ldl -lstdc++ -lgfortran -lmpifort -lopenblas

INCLUDES = -I$(INCPEXSI) -I$(INCSLUDIST) -I$(INCSSPARSE) -I$(INCMUMPS) -I$(INCHYPRE) -I$(INCQHULL) -I$(INCCUDA) -I$(INCMAGMA)

