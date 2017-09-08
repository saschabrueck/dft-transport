# source  /home/seyedb/cp2k-omen/cp2k/cp2k/tools/toolchain/install/setup ; source /home/seyedb/cp2k-omen/omen/cp2komen_envsetup
#
LEX      = flex
YACC     = bison

CXX      = mpic++ 
CC       = gcc 
NVCC     =  

CFLAGS    = -g -w -Wall -fopenmp 
CXXFLAGS  = -std=c++11 $(CFLAGS)
NVCCFLAGS = 

TOP_DIR       = /home/seyedb/cp2k-omen
LIB_TOP       = $(TOP_DIR)/libs
TOOLCHAIN     = /home/seyedb/cp2k-omen/cp2k/cp2k/tools/toolchain

# Common include paths
INCSSPARSE    = -I$(LIB_TOP)/SuiteSparse/include/
INCMUMPS      = -I$(LIB_TOP)/MUMPS/include/
INCHYPRE      = -I$(LIB_TOP)/hypre/include/
INCQHULL      = -I$(LIB_TOP)/qhull/include/libqhull/
INCMAGMA      =
INCSLUDIST    = -I$(TOOLCHAIN)/install/superlu_dist-5.1.2/include/
INCPEXSI      = -I$(TOOLCHAIN)/install/pexsi-0.10.1/include/

INCCUDA       =

# Common library paths
LIBSSPARSE    = -L$(LIB_TOP)/SuiteSparse/static -lumfpack -lamd -lcholmod -lcolamd -lccolamd -lcamd -lsuitesparseconfig 
LIBMUMPS      = -L$(LIB_TOP)/MUMPS/lib -lzmumps -ldmumps -lmumps_common -lpord -lscalapack 
LIBHYPRE      = -L$(LIB_TOP)/hypre/lib -lHYPRE 
LIBQHULL      = -L$(LIB_TOP)/qhull/lib -lqhullstatic 
LIBMAGMA      =
LIBPARMETIS   = -L$(TOOLCHAIN)/install/parmetis-4.0.3/lib -lparmetis -lmetis 
LIBSLUDIST    = -L$(TOOLCHAIN)/install/superlu_dist-5.1.2/lib -lsuperlu_dist 
LIBPEXSI      = -L$(TOOLCHAIN)/install/pexsi-0.10.1/lib -lpexsi 
LIBPARDISO    = -L$(LIB_TOP)/pardiso/lib -lpardiso500-MPI-GNU472-X86-64
LIBCP2K       = -L/home/seyedb/cp2k-omen/cp2k/cp2k/lib/local/popt -lcp2k -lxsmmf -lxsmm -lderiv -lint -lxcf90 -lxc -lfftw3 

LIBCUDA       =

LFLAGS   = $(LIBCP2K) $(LIBPEXSI) $(LIBPARDISO) $(LIBPARMETIS) $(LIBSLUDIST) $(LIBSSPARSE) $(LIBMUMPS) $(LIBHYPRE) $(LIBQHULL) $(LIBCUDA) $(LIBMAGMA)

DFLAGS   = -DAdd_ 

LIBS     = -lrt -ldl -lstdc++ -lgfortran -lmpifort -lopenblas

INCLUDES = $(INCPEXSI) $(INCSLUDIST) $(INCSSPARSE) $(INCMUMPS) $(INCHYPRE) $(INCQHULL) $(INCCUDA) $(INCMAGMA)

