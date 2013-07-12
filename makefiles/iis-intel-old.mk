# Add path to old Intel compiler in /usr/local/mpich2/intel/bin 
# no further environment variables need to be set

# Compiler and MKL library paths
MKLLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t
MKLINC       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/include/em64t
INTLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/lib
INTINC       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/include

# Linear Solvers and Inverters
UMFLIB       =/home/mauro/sw/suitesparse/intel-2011/lib
UMFINC       =/home/mauro/sw/suitesparse/intel-2011/include

# CP2K library path
CP2KLIB = /usr/zupo/home/bruecks/cp2k/cp2k/lib/Intel-Old/popt

CPP = mpicxx
CC =

CFLAGS = -g -wd981 # -wd981 is recommended for c++ even by intel
CXXFLAGS = -std=c++0x -openmp $(CFLAGS)

PARDISO_SO = /home/mauro/sw/pardiso/lib/libpardiso490-INTEL120-X86-64.so

LFLAGS = -L$(CP2KLIB) -L$(INTLIB) -L$(MKLLIB) -L$(UMFLIB)

DFLAGS = -DAdd_ -DMPICH_IGNORE_CXX_SEEK

INCLUDES = -I$(INTINC) -I$(MKLINC) -I$(UMFINC)

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	-lifcore -lifport -lpthread -limf -lsvml -lintlc -lpthread -lm \
	-lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig -lcamd -lccolamd -lmetis \
	$(MKLLIB)/libmkl_scalapack_lp64.a $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a $(MKLLIB)/libmkl_blacs_intelmpi_lp64.a
