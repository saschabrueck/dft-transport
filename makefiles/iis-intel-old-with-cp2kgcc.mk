# Old Intel compiler and CP2K compiled with GNU

LEX           = flex
YACC          = bison

# Compiler and MKL library paths
MKLLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t
MKLINC       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/include/em64t
INTLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/lib
INTINC       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/include

GFORTL       = /home/mauro/sw/gcc/4.8.1/lib64

# Linear Solvers and Inverters
UMFLIB       =/home/mauro/sw/suitesparse/intel-2011/lib
UMFINC       =/home/mauro/sw/suitesparse/intel-2011/include

#METISLIB      = $(LIB_TOP)/Metis/libmetis.a
MUMPSLIB      = $(LIB_TOP)/MUMPS/lib/libzmumps.a
PORDLIB       = $(LIB_TOP)/MUMPS/lib/libpord.a
INCMPS        = -I$(LIB_TOP)/MUMPS/include/
INCPOR        = -I$(LIB_TOP)/MUMPS/PORD/include/

ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a

LIB_TOP       = /usr/ela/home/bruecks/GreenSolver
INCQHU        = -I$(LIB_TOP)/QHULL/src/
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a

INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a

# CP2K library path
#CP2KLIB = /usr/ela/home/bruecks/cp2k/cp2k/lib/iis_ee-x86-64-gcc-openblas-scalapack/popt
CP2KLIB = /usr/ela/home/bruecks/cp2k/cp2k/lib/iis_ee-x86-64-gcc-openblas-scalapack-no_native/popt
#CP2KLIB = /usr/ela/home/bruecks/cp2k2/trunk/cp2k/lib/iis_ee-x86-64-gcc-openblas-scalapack-no_native/popt

CPP = /usr/local/mpich2/intel/bin/mpicxx
#CPP = /usr/local/mpich2-1.5/intel/bin/mpicxx
CC =
GCC = icc

CFLAGS = -g -wd981 # -wd981 is recommended for c++ even by intel
CXXFLAGS = -std=c++0x -openmp $(CFLAGS)

PARDISO_SO = /home/mauro/sw/pardiso/lib/libpardiso490-INTEL120-X86-64.so

LFLAGS = -L$(CP2KLIB) -L$(INTLIB) -L$(MKLLIB) -L$(UMFLIB) -L/usr/local/mpich2/intel/lib/

DFLAGS = -DAdd_ -DMPICH_IGNORE_CXX_SEEK

INCLUDES = -I$(INTINC) -I$(MKLINC) -I$(UMFINC) $(INCQHU) $(INCMPS) $(INCPOR) $(INCAZT)

LIBS = $(QHULLLIB) -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib -Wl,-rpath,$(GFORTL) -L$(GFORTL) -lgfortran \
	$(MUMPSLIB) $(PORDLIB) $(AZTECLIB) $(ARPACKLIB)\
	-Wl,-rpath,$(INTLIB) -lifcore -lifport -lpthread -limf -lsvml -lintlc -lpthread -lm \
	-lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig -lcamd -lccolamd -lmetis \
        -Wl,-rpath,$(MKLLIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64
#	$(MKLLIB)/libmkl_scalapack_lp64.a $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a
