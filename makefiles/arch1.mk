LEX      = flex
YACC     = bison

CPP      = mpic++ 
GCC      = gcc 

CFLAGS   = -g -w -Wall -fopenmp 
CXXFLAGS = -std=c++11 $(CFLAGS)

TOP_DIR        = /data/seyedb
TOOLCHAIN_DIR = /data/vjoost/toolchain-gcc492/install

# Common include paths
IAZTEC       = $(TOP_DIR)/Aztec2.1.1.0/lib/
IPORD        = $(TOP_DIR)/MUMPS_4.10.0/PORD/include/
IMUMPS       = $(TOP_DIR)/MUMPS_4.10.0/include/
ISPARSE      = $(TOP_DIR)/SuiteSparse_config/include/
ISUITESPARSE = $(TOP_DIR)/SuiteSparse_config/
IQHULL       = $(TOP_DIR)/qhull-2012.1/src/libqhull/
ISLU         = $(TOP_DIR)/SuperLU_DIST_3.3/SRC

# Common library paths
LAZTEC       = $(TOP_DIR)/Aztec2.1.1.0/lib
LSLUD        = $(TOP_DIR)/SuperLU_DIST_3.3/lib
LMETIS       = $(TOP_DIR)/parmetis-4.0.2/build/Linux-x86_64/libmetis 
LPARMETIS    = $(TOP_DIR)/parmetis-4.0.2/build/Linux-x86_64/libparmetis
LESMUMPS     = $(TOP_DIR)/scotch_6.0.0_esmumps/lib
LMUMPS       = $(TOP_DIR)/MUMPS_4.10.0/lib 
LSPARSE      = $(TOP_DIR)/SuiteSparse_config/lib
LSUITESPARSE = $(TOP_DIR)/SuiteSparse_config 
LQHULL       = $(TOP_DIR)/qhull-2012.1/lib 
LARPACK      = $(TOP_DIR)/ARPACK
LIBPEXSI     = $(TOP_DIR)/pexsi_v0.5.5/src
#LIBPEXSI     = $(TOP_DIR)/pexsi/src
LIBINT       = $(TOOLCHAIN_DIR)/lib/ 
LIBXC        = $(TOOLCHAIN_DIR)/lib/
SCALAPACK    = $(TOOLCHAIN_DIR)/lib/ 

# CP2K library path
CP2KLIB  = $(TOP_DIR)/cp2k-omen/cp2k/cp2k/lib/local/pdbg/

LFLAGS   = -L$(CP2KLIB) -L$(SCALAPACK) -L$(LAZTEC) -L$(LMETIS) -L$(LPARMETIS) -L$(LSLUD) -L$(LESMUMPS) \
           -L$(LMUMPS) -L$(LSPARSE) -L$(LSUITESPARSE) -L$(LQHULL) -L$(LARPACK) \
           -L$(LIBINT) -L$(LIBXC) -L$(LIBPEXSI)

DFLAGS   = -DAdd_

PARDISO_SO = -L$(TOP_DIR)/pardiso/lib -lpardiso500-MPI-GNU472-X86-64 

LIBS     = -lm -lgfortran -lstdc++ \
           -lcp2k \
           -lfftw3 -laztec -lpexsi \
           -lzmumps -ldmumps -lmumps_common -lpord \
           -lumfpack -lamd -lccolamd -lcholmod -lcolamd -lcamd -lccolamd \
           -lparmetis -lmetis -lesmumps -lsuperlu_dist_3.3 -lscotch -lscotcherr \
           -lscalapack -lreflapack -lrefblas \
           -lsuitesparseconfig \
           -lqhullstatic -larpack -lderiv -lint -lxc -lmpifort
       
INCLUDES = -I$(IAZTEC) -I$(IPORD) -I$(IMUMPS) -I$(ISPARSE) -I$(ISUITESPARSE) -I$(IQHULL) -I$(ISLU)

