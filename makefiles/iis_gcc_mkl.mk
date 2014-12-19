LEX           = flex
YACC          = bison

GFORTL       = /home/nanotcad/sw/gcc/4.8.2/dahu/lib64

LIB_TOP       = /usr/ela/home/bruecks/GreenSolver_Intel

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
CP2KLIB = /home/nanotcad/sw/cp2k/lib/iis_gcc_mkl/popt

CPP = /home/nanotcad/sw/mpich/3.1-gcc-4.8.2/dahu/bin/mpicxx
GCC = icc

CFLAGS = -g -fopenmp
CXXFLAGS = $(CFLAGS)

DMALLOC = -L/usr/local/allinea/tools/lib/64/  -ldmallocthcxx -z muldefs

LFLAGS = -L$(CP2KLIB) -L$(MKLLIB) -L$(INTLIB)

DFLAGS = -DAdd_ -DMPICH_IGNORE_CXX_SEEK -DHAVE_UMFPACK -Dlibcp2k

LINLIN = /home/bruecks/CSelInv/EXAMPLES/C2Finterface.o /home/bruecks/CSelInv/LIB/libcsupldlt.a

INCLUDES = -I$(MKLINC) $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU)

LIBS = -lcp2k \
	$(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(MUMPDLIB) $(MUMPSCOM) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) \
        -Wl,-rpath,$(MKLLIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 \
	-Wl,-rpath,$(INTLIB) -lifcore -lifport -limf -lsvml -lintlc -liomp5 -lpthread \
	-lm -lrt -fopenmp -Wl,-rpath,$(GFORTL) -L$(GFORTL) -lgfortran -lstdc++
