ifneq (,$(findstring .mk,mauro-bh2-intelmkl))
ARCH = ./makeinc.d/local
else
ARCH = ./makeinc.d/local
endif
include $(ARCH)

LDFLAGS = $(LFLAGS) $(DFLAGS)

MYEXECUTABLE = mymain
SOURCES = Utilities.C CSR.C Spike.C Pardiso.C

INTEL_SOURCES = intel_Spike.C
INTEL_FLAGS = -I../libs/intel_spike/include -L../libs/intel_spike/examples/examples_c/../../lib/em64t  -lmkl_solver  -lmkl_lapack -lmkl -lguide -lpthread -lm -lspike_adapt -lspike_adapt_de -lspike_adapt_grid_f -lspike_mpi_comm_intelmpi -lspike 

OBJECTS = $(SOURCES:.C=.o)

all: $(SOURCES) $(MYEXECUTABLE)

.C.o:
	$(CPP) $(INCLUDES) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

$(MYEXECUTABLE): $(MYEXECUTABLE).cc $(SOURCES_C) $(OBJECTS)
	$(CPP) $(CXXFLAGS) $(INCLUDES) $(MYEXECUTABLE).cc -o $@ $(OBJECTS)  $(LDFLAGS) $(LIBS)

intel: $(SOURCES)
	$(CPP) $(CXXFLAGS) $(INCLUDES) intel_mymain.cc intel_Spike.C -o intel_mymain $(LDFLAGS) $(LIBS) $(INTEL_FLAGS) 

.PHONY: clean
clean:
	rm -rf *.o $(MYEXECUTABLE)
