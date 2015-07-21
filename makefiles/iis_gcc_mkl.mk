LEX           = flex
YACC          = bison

GFORTL       = /home/nanotcad/sw/gcc/4.8.2/dahu/lib64

LIB_TOP      = /home/bruecks/GreenSolver_Intel

LIB_SW       = /home/nanotcad/sw

# MKL library paths
MKLLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/lib/em64t
MKLINC       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/mkl/include/em64t
INTLIB       = /usr/pack/intel_compiler-11.1.075-af/Linux-x86_64/lib

# Common include paths
INCAZT        = -I$(LIB_TOP)/AZTEC/lib/
INCQHU        = -I$(LIB_TOP)/QHULL/src/
INCAMD        = -I$(LIB_TOP)/AMD/Include/
INCPDI        = -I$(LIB_TOP)/PDIV/
INCUFC        = -I$(LIB_TOP)/UFconfig/
INCUMF        = -I$(LIB_TOP)/UMFPACK/Include
INCSLU        = -I$(LIB_SW)/SuperLU_DIST_3.3/SRC/
INCPEX        = -I$(LIB_SW)/pexsi_v0.8.0/include/
INCMPS        = -I$(LIB_SW)/MUMPS_5.0.0/include/
INCPOR        = -I$(LIB_SW)/MUMPS_5.0.0/PORD/include/

# Common library paths
ARPACKLIB     = $(LIB_TOP)/ARPACK2/libarpack.a
AMDLIB        = $(LIB_TOP)/AMD/Lib/libamd.a
AZTECLIB      = $(LIB_TOP)/AZTEC/lib/libaztec.a
PDIVLIB       = $(LIB_TOP)/PDIV/libpdiv.a
QHULLLIB      = $(LIB_TOP)/QHULL/src/libqhull.a
UMFPACKLIB    = $(LIB_TOP)/UMFPACK/Lib/libumfpack.a
SLULIB        = $(LIB_SW)/SuperLU_DIST_3.3/lib/libsuperlu_dist_3.3.a
PEXLIB        = $(LIB_SW)/pexsi_v0.8.0/src/libpexsi_dahu.a
PARMLIB       = $(LIB_SW)/parmetis/lib/libparmetis.a $(LIB_SW)/parmetis/lib/libmetis.a
MUMPSLIB      = $(LIB_SW)/MUMPS_5.0.0/lib/libzmumpsdahu.a
MUMPSCOM      = $(LIB_SW)/MUMPS_5.0.0/lib/libmumps_commondahu.a
PORDLIB       = $(LIB_SW)/MUMPS_5.0.0/PORD/lib/libpord.a

# CP2K library path
CP2KLIB = /home/nanotcad/sw/cp2k/lib/iis_gcc_mkl/popt

CPP = /home/nanotcad/sw/mpich/3.1-gcc-4.8.2/dahu/bin/mpicxx
GCC = icc

CFLAGS = -g -fopenmp
CXXFLAGS = $(CFLAGS)

DMALLOC = -L/usr/local/allinea/tools/lib/64/  -ldmallocthcxx -z muldefs

LFLAGS = -L$(CP2KLIB) -L$(MKLLIB) -L$(INTLIB)

DFLAGS = -DAdd_ -DMPICH_IGNORE_CXX_SEEK -DHAVE_MUMPS -Dlibcp2k

INCLUDES = -I$(MKLINC) $(INCSLU) $(INCUFC) $(INCUMF) $(INCMPS) $(INCPOR) $(INCAMD) $(INCAZT) $(INCQHU) $(INCPEX)

LIBS = -lcp2k \
	$(UMFPACKLIB) $(AMDLIB) $(MUMPSLIB) $(MUMPDLIB) $(MUMPSCOM) $(PORDLIB) $(METISLIB) $(AZTECLIB) $(QHULLLIB) $(ARPACKLIB) \
	$(PEXLIB) $(SLULIB) $(PARMLIB) \
        -Wl,-rpath,$(MKLLIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 \
	-Wl,-rpath,$(INTLIB) -lifcore -lifport -limf -lsvml -lintlc -liomp5 -lpthread \
	-lm -lrt -fopenmp -Wl,-rpath,$(GFORTL) -L$(GFORTL) -lgfortran -lstdc++
