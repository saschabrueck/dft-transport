
ARCH = ../makefiles/arch1.mk 
include $(ARCH)

LDFLAGS = $(LFLAGS) $(DFLAGS)

SOURCES = main.cpp cDBCSR_tools.cpp c_scf.cpp 
CFILES = GreenSolver.C array_tools.C Utilities.C CSR.C 

OBJECTS = $(SOURCES:.cpp=.o) $(CFILES) $(TIMINGS)
EXECUTABLE = interface

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

.cpp.o:
	$(CPP) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

#.C.o:
#	$(CPP) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf *o $(EXECUTABLE) output.out csr.*
