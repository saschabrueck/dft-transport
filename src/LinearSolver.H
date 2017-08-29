/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __LINEARSOLVER
#define __LINEARSOLVER

#include "Types.H"
#include "CSR.H"
#include "mpi.h"

template <class T>
class LinearSolver {

public:

    LinearSolver(TCSR<T>* mat,MPI_Comm) {}
    LinearSolver(TCSR<T>* mat,MPI_Comm,int) {}
    LinearSolver(TCSR<T>* mat,int,int,MPI_Comm,MPI_Comm) {}
    LinearSolver(TCSR<T>* mat,double*,int,int,MPI_Comm,MPI_Comm) {}
    LinearSolver() {}

    virtual ~LinearSolver(){}
    virtual void prepare() = 0;
    virtual void prepare(int*,int*,int,int,int*,int) = 0;
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int*,int,int,int*,int) = 0;
    virtual void solve_equation(T* res, T* rhs, int no_rhs) = 0;
	
};

#endif

