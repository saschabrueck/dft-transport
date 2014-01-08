#ifndef INTEL_SPIKE_H_
#define INTEL_SPIKE_H_

#include <vector>
#include "CSR.H"

#include <mpi.h>
#include "spike.h"
#include "intel_spike_help.h"


class intel_SPIKE {
  
public:

  intel_SPIKE(TCSR<double>* matrix, double *RHS, int RHS_col, MPI_Comm communicator,
        std::vector<int> partition_lines, int threads_per_rank);
  void solve_full();
  double *solution;

private:

  // our stuff
  TCSR<double> matrix;
  double *RHS;
  MPI_Comm communicator;
  std::vector<int> partition_lines;
  int threads_per_rank;
  int global_rows;
  int fortran_index;
  int partition_height;
  int num_partitions;


  // intel stuff
  int rank, code, nb_procs, i, j, k;
  array_dim2 f;     // what we call RHS
  int f_col;        // columns in the RHS
  int info;
  void *space;
  spike_param_c_interface pspike;
  matrix_data_c_interface mat;

};

#endif
