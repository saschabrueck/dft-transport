/*!\class GFSolver GFSolver.H
\brief Parent class of all the ballistic NEGF solvers (PDIV and RGF). It is called in GreenFunction.C
*/

#ifndef __GFSOLVER
#define __GFSOLVER

#include <mpi.h>
#include "Types.H"
#include "CSR.H"

class GFSolver {

public:

    GFSolver(TCSR<CPX>*,int,MPI_Comm) {}
    GFSolver() {}

    virtual ~GFSolver(){}
    virtual void solve_equation(CPX*,CPX*,CPX*,CPX*,CPX*,CPX*,int*,int*,int,int,\
				int,int,int,int,int*,int*,int) = 0;
	
};

#endif
