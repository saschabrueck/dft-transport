LEX           = flex
YACC          = bison

GFORTL       = /home/mauro/sw/gcc/4.8.1/lib64

LIB_TOP       = /usr/ela/home/bruecks/GreenSolver

# MKL library paths
MKLLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t
MKLINC       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/include/em64t
INTLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/lib

# Common include paths
INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/QHULL/src/
INCAMD        = -I$(LIB_TOP)/AMD/Include/
INCSLU        = -I$(LIB_TOP)/SuperLU_DIST_2.0/SRC/
INCMPS        = -I$(LIB_TOP)/MUMPS_4.10.0/include/
INCPDI        = -I$(LIB_TOP)/PDIV/
INCPOR        = -I$(LIB_TOP)/MUMPS_4.10.0/PORD/include/
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
METISLIB      = $(LIB_TOP)/Metis/libmetis.a
MUMPSLIB      = $(LIB_TOP)/MUMPS_4.10.0/lib/libzmumps.a
MUMPDLIB      = $(LIB_TOP)/MUMPS_4.10.0/lib/libdmumps.a
MUMPSCOM      = $(LIB_TOP)/MUMPS_4.10.0/lib/libmumps_common.a
PDIVLIB       = $(LIB_TOP)/PDIV/libpdiv.a
PORDLIB       = $(LIB_TOP)/MUMPS_4.10.0/lib/libpord.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a

# CP2K library path
CP2KLIB = /usr/ela/home/bruecks/cp2k1/lib/iis_gcc_mkl/psmp

CPP = /usr/local/mpich2-1.5/gcc/bin/mpicxx
GCC = icc

CFLAGS = -g -fopenmp
CXXFLAGS = $(CFLAGS)

LFLAGS = -L$(CP2KLIB) -L$(MKLLIB) -L$(INTLIB)

DFLAGS = -DAdd_ -DMPICH_IGNORE_CXX_SEEK

PARDISO_SO = -L/usr/ela/home/bruecks/Downloads -lpardiso500-MPI-GNU472-X86-64 -Wl,-rpath=/usr/ela/home/bruecks/Downloads

LINLIN = /home/bruecks/CSelInv/EXAMPLES/C2Finterface.o /home/bruecks/CSelInv/LIB/libcsupldlt.a

INCLUDES = -I$(MKLINC) $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU)

LIBS = -lcp2k_lib -lcp2k_base_lib -lcp2k_dbcsr_lib -lcp2k_fft_lib -lcp2k_ma_lib -lcp2k_elpa_lib \
	$(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(MUMPDLIB) $(MUMPSCOM) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) \
        -Wl,-rpath,$(MKLLIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 \
	-Wl,-rpath,$(INTLIB) -lifcore -lifport -limf -lsvml -lintlc -liomp5 -lpthread \
	-lm -lrt -fopenmp -Wl,-rpath,$(GFORTL) -L$(GFORTL) -lgfortran -lstdc++
