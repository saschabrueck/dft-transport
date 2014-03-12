#ifndef __SPIKE_SOLVER_H
#define __SPIKE_SOLVER_H

#include <vector>
#include <mpi.h>
#include "LinearSolver.H"
#include "Spike.H"

#include <iostream> // DEBUG

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
};

template <typename T>
void SpikeSolver<T>::solve_equation(T* result, T* rhs, int no_rhs) {
  Spike<T>* _spike = new Spike<T>(_matrix, rhs, no_rhs, _communicator);
  _spike->solve_full(result);
};
  

#endif
