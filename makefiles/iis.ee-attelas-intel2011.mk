# Assumes you use the stuff under /home/mauro/sw
#
# To set up the environment accordingly:
#
# 	bash
# 	source ~mauro/sw/set_environment-intel2011-openmpi.sh
#
# Or perform the required steps manually for your csh.

MBASE = /home/mauro/sw
IBASE = /usr/pack/intel_compiler-11.1.075-af

CPP = mpic++
CC =

#CFLAGS = -g -Wall -DAdd_ -wd981			# -wd981 is recommended for c++ even by intel
CFLAGS = -g -DAdd_ -wd981 # -wd981 is recommended for c++ even by intel
CXXFLAGS = -std=c++0x -openmp $(CFLAGS)
PARDISO_SO = /home/mauro/sw/pardiso/lib/libpardiso490-INTEL120-X86-64.so
INTEL_LIB = -L$(IBASE)/Linux-x86_64/mkl/lib/em64t \
							-lmkl_intel_lp64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
							-lmkl_sequential -lmkl_core \
							-Wl,-rpath,$(IBASE)/Linux-x86_64/mkl/lib/em64t \
						-L$(IBASE)/Linux-x86_64/lib -lifcore -lifport -lpthread -limf \
							-lsvml -lintlc -lm -lpthread -lm \
							-Wl,-rpath,$(IBASE)/Linux-x86_64/lib

CP2KLIB = /home/mauro/sw/cp2k/cp2k-transport/lib/iis_ee-x86-64-oldintel/popt

# Use scalapack, suitesparse/umfpack is available on this system.
LFLAGS = -L$(CP2KLIB) -L$(OLD_INTEL) -L$(MBASE)/gcc/4.8.1/lib64 -L$(MBASE)/fftw/intel-2011/lib -L$(MBASE)/suitesparse/intel-2011/lib

DFLAGS =

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	 $(INTEL_LIB) \
   -lm \
   -lfftw3 \
   -I$(MBASE)/suitesparse/intel-2011/include -lumfpack -lamd \
	 		-lcholmod -lcolamd -lrt -lsuitesparseconfig \
			-lcamd -lccolamd -lmetis \
	 -lmpi_f77 -Wl,-rpath,$(MBASE)/openmpi/intel-2011/lib

PARDISO_LIBS = $(PARDISO_INTEL_LIB) $(PARDISO_SO) -limf -lsvml -liomp5 -lintlc
