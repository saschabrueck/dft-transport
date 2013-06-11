CPP = mpic++
CC =

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 $(CFLAGS) -DMMC

CP2KLIB = -L/home/mauro/omendft/external/cp2k/cp2k/lib/mauro-bh2-gfortran-scalapack/popt

# Use scalapack, suitesparse/umfpack is available on this system.
LFLAGS = $(CP2KLIB)

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
   -lm -lgfortran -lstdc++ \
   -lfftw3 \
   -lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig \
	 -lscalapack
