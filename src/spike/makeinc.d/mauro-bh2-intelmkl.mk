CPP = mpic++
CC =

CFLAGS = -g -Wall -DAdd_
CXXFLAGS = -std=c++11 $(CFLAGS)

# MKL includes blas and scalapack, suitesparse/umfpack is available on this
# system. OpenMPI needs to be linked explicitly, mpicxx doesn't handle it.
LFLAGS = $(CP2KLIB) \
	 -L/opt/intel/mkl/lib/intel64

DFLAGS = 

INCLUDES = -I/opt/intel/composerxe/mkl/include

LIBS = -lmkl_gf_lp64 -lmkl_scalapack_lp64 -lmkl_sequential -lmkl_core \
       -lmkl_blacs_openmpi_lp64 -lmpi_f77 -lmkl_gnu_thread 
