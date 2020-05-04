# source  /home/hossein/Projects/cp2k-omen/cp2k/tools/toolchain/install/setup ; source /home/hossein/Projects/cp2k-omen/omen/install/cp2komen_envsetup
#
LEX      = flex
YACC     = bison

CXX      = mpic++ 
CC       = gcc 
NVCC     = 

CFLAGS    = -g -w -Wall -fopenmp 
CXXFLAGS  = -std=c++11 $(CFLAGS)
NVCCFLAGS = 

TOP_DIR       = /home/hossein/Projects/cp2k-omen/omen/install
LIB_TOP       = $(TOP_DIR)/libs
TOOLCHAIN     = /home/hossein/Projects/cp2k-omen/cp2k/tools/toolchain

# Common include paths
INCSSPARSE    = -I$(LIB_TOP)/SuiteSparse/include/
INCMUMPS      = -I$(LIB_TOP)/MUMPS/include/
INCHYPRE      = -I$(LIB_TOP)/hypre/include/
INCQHULL      = -I$(LIB_TOP)/qhull/include/libqhull/
INCMAGMA      = 
INCSLUDIST    = -I$(TOOLCHAIN)/install/superlu_dist-6.1.0/include/
INCPEXSI      = -I$(TOOLCHAIN)/install/pexsi-1.2.0/include/

INCCUDA       = 

# Common library paths
LIBSSPARSE    = -L$(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig
LIBMUMPS      = -L$(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack
LIBHYPRE      = -L$(LIB_TOP)/hypre/lib -lHYPRE
LIBQHULL      = -L$(LIB_TOP)/qhull/lib -lqhullstatic
LIBMAGMA      = 
LIBPARMETIS   = -L$(LIB_TOP)/parmetis-4.0.3/lib -lparmetis -lmetis
LIBSLUDIST    = -L$(TOOLCHAIN)/install/superlu_dist-6.1.0/lib -lsuperlu_dist
LIBPEXSI      = -L$(TOOLCHAIN)/install/pexsi-1.2.0/lib -lpexsi
LIBPARDISO    = 
LIBDBCSR      = -L/home/hossein/Projects/cp2k-omen/cp2k/lib/local/popt/exts/dbcsr -ldbcsr
LIBCP2K       = -L/home/hossein/Projects/cp2k-omen/cp2k/lib/local/popt -lcp2k -lxsmmf -lxsmm -lxcf03 -lxc -lint2 -lsymspg -lelpa -lfftw3

LIBCUDA       = 

LFLAGS   = $(LIBCP2K) $(LIBDBCSR) $(LIBPEXSI) $(LIBPARDISO) $(LIBPARMETIS) $(LIBSLUDIST) $(LIBSSPARSE) $(LIBMUMPS) $(LIBHYPRE) $(LIBQHULL) $(LIBCUDA) $(LIBMAGMA)

DFLAGS   = -DAdd_ 

LIBS     = -lrt -ldl -lstdc++ -lgfortran -lmpifort -lopenblas

INCLUDES = $(INCPEXSI) $(INCSLUDIST) $(INCSSPARSE) $(INCMUMPS) $(INCHYPRE) $(INCQHULL) $(INCCUDA) $(INCMAGMA)

