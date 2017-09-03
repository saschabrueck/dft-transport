# source /home/seyedb/cp2k-omen/cp2k/cp2k/tools/toolchain/install/setup
#
LEX      = flex
YACC     = bison

CPP      = mpic++ 
GCC      = gcc 
NVCC     = 

CFLAGS    = -g -w -Wall -fopenmp 
CXXFLAGS  = -std=c++11 $(CFLAGS)
NVCCFLAGS = 

TOP_DIR       = /home/seyedb/cp2k-omen
LIB_TOP       = $(TOP_DIR)/libs
TOOLCHAIN     = $(TOP_DIR)/cp2k/cp2k/tools/toolchain

# Common include paths
INCMUMPS      = $(LIB_TOP)/MUMPS/include/
INCHYPRE      = $(LIB_TOP)/hypre/include/
INCQHULL      = $(LIB_TOP)/qhull/include/libqhull/
INCSSPARSE    = $(LIB_TOP)/SuiteSparse/include/
INCPEXSI      = $(TOOLCHAIN)/install/pexsi-0.10.1/include/
INCSLUDIST    = $(TOOLCHAIN)/install/superlu_dist-5.1.2/include/

# Common library paths
LIBPARMETIS   = $(TOOLCHAIN)/install/parmetis-4.0.3/lib
LIBMUMPS      = $(LIB_TOP)/MUMPS/lib
LIBHYPRE      = $(LIB_TOP)/hypre/lib
LIBQHULL      = $(LIB_TOP)/qhull/lib
LIBSSPARSE    = $(LIB_TOP)/SuiteSparse/static
LIBSLUDIST    = $(TOOLCHAIN)/install/superlu_dist-5.1.2/lib
LIBPEXSI      = $(TOOLCHAIN)/install/pexsi-0.10.1/lib
LIBCP2K       = $(TOP_DIR)/cp2k/cp2k/lib/local/popt

LFLAGS   = -L$(LIBCP2K) -L$(LIBSSPARSE) -L$(LIBPEXSI) -L$(LIBPARMETIS) -L$(LIBMUMPS)\
	-L$(LIBSSPARSE) -L$(LIBSLUDIST) -L$(LIBHYPRE) -L$(LIBQHULL)

DFLAGS   = -DAdd_ 

LIBS     = -lrt -ldl -lstdc++ -lcp2k -lfftw3 -lpexsi \
	-lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig \
	-lparmetis -lmetis -lzmumps -ldmumps -lmumps_common -lpord -lsuperlu_dist \
	-lHYPRE -lscalapack -llapack -lopenblas -lqhullstatic \
	-lxsmmf -lxsmm -lderiv -lint -lxcf90 -lxc -lgfortran -lmpifort

INCLUDES = -I$(INCSLUDIST) -I$(INCHYPRE) -I$(INCSSPARSE) -I$(INCMUMPS) -I$(INCQHULL) -I$(INCPEXSI)

