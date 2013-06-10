CPP = mpicxx
CC =

CFLAGS = -g -Wall
CXXFLAGS = -std=c++11 $(CFLAGS) -DMMC

CP2KLIB = -L/home/mauro/omendft/external/cp2k/cp2k/lib/mauro-bh2-gfortran/popt

# MKL includes blas and scalapack, suitesparse/umfpack is available on this
# system
LFLAGS = $(CP2KLIB) \
	 -L/opt/intel/mkl/lib/intel64

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
   -lm -lgfortran -lstdc++ \
   -lfftw3 -lderiv \
   -lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig \
   -lmkl_gf_lp64 -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64

# test if this is neccessary
#INCLUDES = -I/opt/intel/mkl/include
