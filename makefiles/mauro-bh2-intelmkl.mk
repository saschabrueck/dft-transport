CPP = mpic++
CC =

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 $(CFLAGS) -DMMC

CP2KLIB = -L/home/mauro/omendft/external/cp2k/cp2k/lib/mauro-bh2-gfortran-intelmkl/popt

# MKL includes blas and scalapack, suitesparse/umfpack is available on this
# system. OpenMPI needs to be linked explicitly, mpicxx doesn't handle it.
LFLAGS = $(CP2KLIB) \
	 -L/opt/intel/mkl/lib/intel64

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
   -lm -lgfortran -lstdc++ \
   -lfftw3 \
   -lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig \
   -lmkl_gf_lp64 -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 \
	 -lmpi_f77
# god knows why -lmpi_f77 is needed but the missing references were in 
