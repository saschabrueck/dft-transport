#SCOREP = scorep

CPP = $(SCOREP) CC
GCC = cc

#CFLAGS = -Wall -fopenmp -O3 -ffast-math -funroll-loops -march=native
CFLAGS = -Wall -fopenmp -g
CXXFLAGS = -std=c++11 $(CFLAGS)

LEX           = flex
YACC          = bison

# Directory of the libraries used by OMEN
LIB_TOP       = /project/s503/omendft_libraries
LIB_TOP2      = /project/s503/OMEN_XC30

# Common include paths
INCAZT        = -I$(LIB_TOP2)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP2)/QHULL/src/
INCAMD        = -I$(LIB_TOP2)/AMD/Include/
INCUFC        = -I$(LIB_TOP2)/UFconfig/
INCUMF        = -I$(LIB_TOP2)/UMFPACK/Include/
INCPEX        = -I$(LIB_TOP)/pexsi_v0.6.0/include/

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AZTECLIB      = $(LIB_TOP2)/AZTEC/lib/libaztec.a
QHULLLIB      = $(LIB_TOP2)/QHULL/src/libqhull.a
AMDLIB        = $(LIB_TOP2)/AMD/Lib/libamd.a
UMFPACKLIB    = $(LIB_TOP2)/UMFPACK/Lib/libumfpack.a
PEXSILIB      = $(LIB_TOP)/pexsi_v0.6.0/src/libpexsi_daint.a

LINLIN = $(LIB_TOP)/CSelInv/EXAMPLES/C2Finterface.o $(LIB_TOP)/CSelInv/LIB/libcsupldlt.a

DMALLOC = -L/apps/rosa/ddt/4.1.1/lib/64/ -ldmallocthcxx -z muldefs

LFLAGS = -L$(LIB_TOP)/cp2k/lib/CRAY-XC30-gfortran/psmp/
DFLAGS = -DAdd_
LIBS = $(UMFPACKLIB) $(AMDLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) $(PEXSILIB) \
	-lzmumps -lsuperlu_dist \
	-lcp2k /project/cray/alazzaro/libsmm/affinity/sandybridge_gcc_4.9.0/lib/libsmm_dnn_cray.gnu.a \
        -lfftw3 -lfftw3_threads \
        $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group \
        $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a \
        $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -ldl \
	-lm -lrt -fopenmp -lgfortran -lstdc++ 
INCLUDES = $(INCUFC) $(INCUMF) $(INCAMD) $(INCAZT) $(INCQHU) $(INCPEX)
