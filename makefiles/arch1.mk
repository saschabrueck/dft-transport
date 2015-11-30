LEX      = flex
YACC     = bison

CPP      = mpic++ 
GCC      = gcc 

CFLAGS   = -g -w -Wall -fopenmp 
CXXFLAGS = -std=c++11 $(CFLAGS)

TOP_DIR       = /data/seyedb/cp2k-omen
TOOLCHAIN_DIR = /data/vjoost/toolchain-r16358/install
LIB_TOP       = $(TOP_DIR)/libs

# Common include paths
INCMUMPS      = $(LIB_TOP)/mumps/include/
INCAZTEC      = $(LIB_TOP)/aztec/include/
INCQHULL      = $(LIB_TOP)/qhull/include/
INCSSPARSE    = $(LIB_TOP)/suitesparse/include/
INCPEXSI      = $(LIB_TOP)/pexsi/include/
INCSLUDIST    = $(TOOLCHAIN_DIR)/include/superlu_dist_3.3/

# Common library paths
LIBARPACK     = $(LIB_TOP)/arpack/lib
LIBMUMPS      = $(LIB_TOP)/mumps/lib
LIBAZTEC      = $(LIB_TOP)/aztec/lib
LIBQHULL      = $(LIB_TOP)/qhull/lib
LIBSSPARSE    = $(LIB_TOP)/suitesparse/lib
LIBCP2K       = $(TOP_DIR)/cp2k/cp2k/lib/local_omencp2k/popt
TOOLCHAIN_LIB = $(TOOLCHAIN_DIR)/lib

LFLAGS   = -L$(LIBCP2K) -L$(TOOLCHAIN_LIB) -L$(LIBAZTEC) \
           -L$(LIBMUMPS) -L$(LIBSSPARSE) -L$(LIBQHULL) -L$(LIBARPACK)

DFLAGS   = -Dlibcp2k -DAdd_

LIBS     = -lm -lgfortran -lstdc++ \
           -lcp2k \
           -lfftw3 -laztec -lpexsi_linux_v0.9.0 \
           -lumfpack -lamd -lccolamd -lcholmod -lcolamd -lcamd -lccolamd \
           -lparmetis -lmetis -lzmumps -ldmumps -lmumps_common -lpord -lsuperlu_dist_3.3 \
           -lscalapack -lreflapack -lrefblas \
           -lsuitesparseconfig \
           -lqhullstatic -larpack -lderiv -lint -lxcf90 -lxc -lmpifort

INCLUDES = -I$(INCAZTEC) -I$(INCSLUDIST) -I$(INCSSPARSE) -I$(INCMUMPS) -I$(INCQHULL) -I$(INCPEXSI)

