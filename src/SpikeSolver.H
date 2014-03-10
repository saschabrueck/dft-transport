#ifndef __SPIKE_SOLVER_H
#define __SPIKE_SOLVER_H

#include <vector>
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
  void prepare(int bandwidth, std::vector<int> partition_lines);
  void solve_equation(T* res, T* rhs, int no_rhs);

 private:
  TCSR<T>*          _matrix;
  MPI_Comm          _communicator;
  int               _bandwidth;
  std::vector<int>  _partition_lines;
};

template <typename T>
void SpikeSolver<T>::prepare(int bandwidth, std::vector<int> partition_lines) {
  _bandwidth = bandwidth;
  _partition_lines = partition_lines;
};

template <typename T>
SpikeSolver<T>::SpikeSolver(TCSR<T>* matrix, MPI_Comm communicator)
                       : _matrix(matrix),
                         _communicator(communicator) {
};

template <typename T>
void SpikeSolver<T>::solve_equation(T* result, T* rhs, int no_rhs) {
  // Change later:
  bool use_dense_for_diagonal = false;
  Spike<T>* _spike = new Spike<T>(_matrix, _bandwidth, rhs, no_rhs, 
                                  _communicator, _partition_lines, 
                                  use_dense_for_diagonal);
  _spike->solve_full(result);
};
  

#endif
