CPP = mpicxx
CC =

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 $(CFLAGS)

CP2KLIB = -L/data/seyedb/clean/cp2k/cp2k/lib/Linux-x86-64-gfortran/pdbg/

LFLAGS = $(CP2KLIB) \
   -L/data/seyedb/scalapack/install/lib/ \
   -L/data/vjoost/libint_ham/install/lib/ \
   -L/data/vjoost/libxc-2.0.1/install/lib/ \
   -L/data/seyedb/SuiteSparse_config/lib #\
   -L/opt/intel/composerxe-2011.3.174/mkl/lib/intel64

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
   -lm -lgfortran -lstdc++ \
   -lfftw3 -lderiv -lint -lxc \
   -lscalapack -lreflapack -lrefblas \
   -lumfpack -lamd -lcholmod -lcolamd -lrt -lsuitesparseconfig #\
   -lmkl_gf_lp64 -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64

INCLUDES = -I/data/seyedb/SuiteSparse_config/include
