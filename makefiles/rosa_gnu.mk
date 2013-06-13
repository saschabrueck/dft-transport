#CPP = g++
#CPP = mpicc
CPP = CC
CC = 

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 $(CFLAGS)

#--------------------------------------------------------------------------------------------
# modify paths to the following libraries according to your local installations :
# CP2K, ScaLapack, BLAS and LAPACK, libint, libxc (use the same libraries linked to CP2K)
# UMFPACK, SuiteSparse_config, AMD, CHOLMOD, COLAMD
#--------------------------------------------------------------------------------------------

# Top library directory, relative to OMEN
LIB_TOP       = /users/bruecks/interface

# Common include paths
INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/QHULL/src/
INCAMD        = -I$(LIB_TOP)/AMD/Include/
INCSLU        = -I$(LIB_TOP)/SuperLU_DIST_2.0/SRC/
INCMPS        = -I$(LIB_TOP)/MUMPS/include/
INCPDI        = -I$(LIB_TOP)/PDIV/
INCPOR        = -I$(LIB_TOP)/MUMPS/PORD/include/
INCSLU        = -I$(LIB_TOP)/SuperLU_DIST_2.0/SRC
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK/libarpack.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
METISLIB      = $(LIB_TOP)/Metis/libmetis.a
MUMPSLIB      = $(LIB_TOP)/MUMPS/lib/libzmumps.a
PDIVLIB       = $(LIB_TOP)/PDIV/libpdiv.a
PORDLIB       = $(LIB_TOP)/MUMPS/lib/libpord.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
SUPERLULIB    = $(LIB_TOP)/SuperLU_DIST_2.0/Lib/libsuperlu_dist_2.0.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a

LFLAGS = -L/users/bruecks/cp2k/lib/CRAY-XE6-gfortran-hwtopo/popt/ #\
	-L/opt/xt-libsci/default/gnu/47/interlagos/lib #\
         -L/users/bruecks/GreenSolver/UMFPACK/Lib/ #\
	-L/data/seyedb/scalapack/install/lib/ #\
	-L/data/vjoost/libint_ham/install/lib/ #\
	-L/data/vjoost/libxc-1.2.0/install/lib/ #\
	-L/data/seyedb/SuiteSparse_config/lib #\
	-L/opt/intel/composerxe-2011.3.174/mkl/lib/intel64
DFLAGS =
LIBS = $(SUPERLULIB) -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lm -lgfortran -lstdc++\
	-lfftw3 -lsci_gnu #\
        -lumfpack #\
	-lscalapack -lreflapack -lrefblas #\
        -lderiv -lint -lxc #\
	-lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig #\
	-lmkl_gf_lp64 -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 
#INCLUDES = -I/data/seyedb/SuiteSparse_config/include 
#INCLUDES = -I/users/bruecks/GreenSolver/UMFPACK/Include/ \
	   -I/users/bruecks/GreenSolver/UFconfig/
#INCLUDES = -I/opt/xt-libsci/default/gnu/47/interlagos/include/
INCLUDES = $(INCSLU) $(INCUMF)
