#ifndef __LINEARSOLVER
#define __LINEARSOLVER

#include <vector>
#include "Types.H"
#include "CSR.H"
#include "mpi.h"

template <class T>
class LinearSolver {

public:

    LinearSolver(TCSR<T>* mat,MPI_Comm) {}
    LinearSolver(TCSR<T>* mat,MPI_Comm,int) {}
    LinearSolver() {}

    virtual ~LinearSolver(){}
    virtual void prepare() = 0;
    virtual void prepare(int*,int*,int,int,int*,int) = 0;
    virtual void solve_equation(T* res, T* rhs, int no_rhs) = 0;
	
};

#endif

