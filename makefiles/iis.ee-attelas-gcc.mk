# assumes you use the stuff under /home/mauro/sw
#
# TODO: ensure that there's no problem for pardiso when linking against a newer
# limf.so than the one it was compiled with (Mathieu's intel)

MBASE = /home/mauro/sw

CPP = mpic++
CC =

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 -fopenmp $(CFLAGS)
#PARDISO_SO = /home/mauro/sw/pardiso/lib/libpardiso490-INTEL120-X86-64.so
PARDISO_SO = /home/mauro/sw/pardiso/lib/libpardiso491-GNU430-X86-64.so
PARDISO_INTEL_LIB = -L/usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t/ -L/usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/lib/ -lifcore -lifport -lpthread -limf -lsvml -lintlc -lm -lpthread -lm -Wl,-rpath,/usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/lib/ -Wl,-rpath,/usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t

CP2KLIB = /usr/zupo/home/bruecks/cp2k/cp2k/lib/iis_ee-x86-64-gcc-openblas-scalapack/popt

# Use scalapack, suitesparse/umfpack is available on this system.
LFLAGS = -L$(CP2KLIB) -L$(MBASE)/gcc/4.8.1/lib64 -L$(MBASE)/fftw/gcc-attelas/lib -L$(MBASE)/suitesparse/gcc-attelas/lib -L$(MBASE)/openblas/gcc-attelas/lib -L$(MBASE)/scalapack/gcc-attelas/lib

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
   -lm -lgfortran -lstdc++ \
   -lfftw3 \
   -I$(MBASE)/suitesparse/gcc-attelas/include -lumfpack -lamd \
	 		-lcholmod -lcolamd -lrt -lsuitesparseconfig \
			-lcamd -lccolamd -lmetis \
	 -lopenblas \
	 -lscalapack

PARDISO_LIBS = $(PARDISO_INTEL_LIB) $(PARDISO_SO) -limf -lsvml -liomp5 -lintlc

LINLIN = /home/bruecks/CSelInv/EXAMPLES/C2Finterface.o /home/bruecks/CSelInv/LIB/libcsupldlt.a

INCLUDES= -I$(MBASE)/suitesparse/gcc-attelas/include
