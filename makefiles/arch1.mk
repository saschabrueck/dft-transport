LEX      = flex
YACC     = bison

CPP      = mpicxx
GCC      = gcc 

CFLAGS   = -g -Wall -fopenmp 
CXXFLAGS = -std=c++11 $(CFLAGS)

TOP_DIR  = /data/seyedb

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
LMETIS       = $(TOP_DIR)/parmetis-4.0.2/build/Linux-x86_64/libmetis 
LPARMETIS    = $(TOP_DIR)/parmetis-4.0.2/build/Linux-x86_64/libparmetis
LESMUMPS     = $(TOP_DIR)/scotch_6.0.0_esmumps/lib
LMUMPS       = $(TOP_DIR)/MUMPS_4.10.0/lib 
LSPARSE      = $(TOP_DIR)/SuiteSparse_config/lib
LSUITESPARSE = $(TOP_DIR)/SuiteSparse_config 
LQHULL       = $(TOP_DIR)/qhull-2012.1/lib 
LARPACK      = $(TOP_DIR)/ARPACK
LIBINT       = /data/vjoost/libint_ham/install/lib/ 
LIBXC        = /data/vjoost/libxc-2.0.1/install/lib/

SCALAPACK = $(TOP_DIR)/scalapack/install/lib/ 

# CP2K library path
CP2KLIB  = $(TOP_DIR)/cp2k/cp2k/lib/Linux-x86-64-gfortran/pdbg/

LFLAGS   = -L$(CP2KLIB) -L$(SCALAPACK) -L$(LAZTEC) -L$(LMETIS) -L$(LPARMETIS) -L$(LESMUMPS) \
           -L$(LMUMPS) -L$(LSPARSE) -L$(LSUITESPARSE) -L$(LQHULL) -L$(LARPACK) -L$(LIBINT) -L$(LIBXC)

DFLAGS   = -DAdd_

PARDISO_SO = -L$(TOP_DIR)/pardiso/lib -lpardiso500-MPI-GNU472-X86-64 

LIBS     = -lm -lgfortran -lstdc++ \
           -lcp2kstart -lcp2k -lcp2kinput -lcp2kbase -lcp2kmpiwrap -lcp2kpw -lcp2kdbcsrwrap -ldbcsr -lcp2kacc \
           -lcp2kfm -lcp2kcommon -lcp2ktmc -lcp2kao -lcp2kxc -lcp2kma -lcp2kfft -lcp2kgrid -lcp2kmetadyn_tools -lcp2kmachine \
           -lfftw3 -laztec \
           -lzmumps -ldmumps -lmumps_common -lpord \
           -lumfpack -lamd -lccolamd -lcholmod -lcolamd -lcamd -lccolamd \
           -lparmetis -lmetis -lesmumps -lscotch -lscotcherr \
           -lscalapack -lreflapack -lrefblas \
           -lsuitesparseconfig \
           -lqhull -larpack -lderiv -lint -lxc
       
INCLUDES = -I$(IAZTEC) -I$(IPORD) -I$(IMUMPS) -I$(ISPARSE) -I$(ISUITESPARSE) -I$(IQHULL) -I$(ISLU)

