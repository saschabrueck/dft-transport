#ifndef __INJECTION
#define __INJECTION

#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Types.H"
#include "Utilities.H"
#include "Blas.H"
#include "ScaLapack.H"
#include "Connection.H"
#include "LinearSolver.H"
#include "Umfpack.H"

class Injection{
public:
    virtual int execute(TCSR<CPX>**,int,int,CPX*,CPX*,int,MPI_Comm){};
};

#endif
