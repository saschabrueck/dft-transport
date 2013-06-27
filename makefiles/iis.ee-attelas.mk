# assumes you use the stuff under /home/mauro/sw

MBASE = /home/mauro/sw

CPP = mpic++
CC =

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 $(CFLAGS) -DMMC

CP2KLIB = /home/mauro/sw/cp2k/cp2k/lib/iis_ee-x86-64-gcc-openblas-scalapack/popt

# Use scalapack, suitesparse/umfpack is available on this system.
LFLAGS = -L$(CP2KLIB) -L$(MBASE)/gcc/4.8.1/lib64 -L$(MBASE)/fftw/gcc/lib -L$(MBASE)/suitesparse/lib -L$(MBASE)/openblas/lib -L$(MBASE)/scalapack/lib

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
   -lm -lgfortran -lstdc++ \
   -lfftw3 \
   -I$(MBASE)/suitesparse/include -lumfpack -lamd \
	 		-lcholmod -lcolamd -lrt -lsuitesparseconfig \
			-lcamd -lccolamd -lmetis \
	 -lopenblas \
	 -lscalapack
