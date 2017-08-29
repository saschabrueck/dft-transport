/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __SPIKE_SOLVER_H
#define __SPIKE_SOLVER_H

#include <vector>
#include <mpi.h>
#include "LinearSolver.H"
#include "Spike.H"

#include <omp.h>

template <typename T>
class SpikeSolver : public LinearSolver<T> {
 public:
  SpikeSolver(TCSR<T>* matrix, MPI_Comm communicator);
  ~SpikeSolver() {};

  // the LinearSolver interface
  void prepare() {};
  void prepare(int*, int*, int, int, int*, int) {};
  void solve_equation(T* res, T* rhs, int no_rhs);

 private:
  TCSR<T>*          _matrix;
  MPI_Comm          _communicator;
};

template <typename T>
SpikeSolver<T>::SpikeSolver(TCSR<T>* matrix, MPI_Comm communicator)
                       : _matrix(matrix),
                         _communicator(communicator) {
  // DEBUG
  printf("SpikeSolver: using %d threads\n", omp_get_max_threads());
};

template <typename T>
void SpikeSolver<T>::solve_equation(T* result, T* rhs, int no_rhs) {
  Spike<T>* _spike = new Spike<T>(_matrix, rhs, no_rhs, _communicator);
  _spike->solve_full(result);
};
  

#endif
