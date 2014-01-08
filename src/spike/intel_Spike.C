#include "intel_Spike.H"

/** \brief Constructor
 *
 *  Initializes the system on each MPI rank.
 *
 *  \param[in] matrix   The horizontal partition of the matrix describing the
 *         system of linear equations in CSR format for the respective rank.
 *         The partitions are assumed to be distributed according to the order
 *         of the rank, that is the first partition is held by rank 0 in the
 *         communicator, the second by rank 1 and so on.
 *
 *  \param[in] RHS      The horizontal partition of the right hand side as 
 *         array of size row * RHS_col where the number of rows is determined
 *         by the partition_lines argument (see below). The NULL pointer 
 *         indicates an all zero RHS for this partition or a RHS that is yet to
 *         be computed.
 *
 *  \param[in] RHS_col  The number of right hand sides (i.e. the number of
 *         columns for the right hand side). Rows will be extracted from
 *         partition_lines.
 *
 *  \param[in] communicator  The MPI communicator to be used for solving the
 *         system.
 *
 *  \param[in] partition_lines  A vector containing the row number of the first
 *         line of the each horizontal partition of the systems, starting with 0
 *         and ending with the total number of rows + 1.
 *
 */
intel_SPIKE::intel_SPIKE(TCSR<double> *matrix, double *RHS, int RHS_col,
             MPI_Comm communicator, std::vector<int> partition_lines,
             int threads_per_rank):
                  matrix(matrix), RHS(RHS), communicator(communicator),
                  partition_lines(partition_lines), 
                  threads_per_rank(threads_per_rank), f_col(RHS_col) 
{

  MPI_Comm_size(communicator, &nb_procs);
  MPI_Comm_rank(communicator, &rank);
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  solution = NULL;
  num_partitions = partition_lines.size()-1;
  global_rows = partition_lines[num_partitions];
  fortran_index = matrix->findx;
  
  partition_height = partition_lines[rank+1]-partition_lines[rank];

  // intel SPIKE stuff
  pspike.nbprocs = nb_procs;
  pspike.rank = rank;
  spike_default(&pspike);
  pspike.tp = 1;          // we have it already distributed
  pspike.autoadapt = 0;   // auto adapt = no
  pspike.timing = 1;      // show some timing
  pspike.comd = 0;        // don't be verbose

  // global information
  mat.format = 'S';       // set to CSR here
  mat.ASTRU = 'G';
  mat.DIAGDO = 'N';
  mat.n = global_rows;    // how large the matrix is globally

  // intel SPIKE equivalent of our partition_lines:
  alloc1D_I(&(mat.sizeA), num_partitions);
  space = mat.sizeA.base;
  for (int i = 0; i < num_partitions; ++i) {
    int part = i + 1;
    int part_size = partition_lines[i + 1] - partition_lines[i];
    setElem1D_I(&(mat.sizeA), part, part_size);
  }

  // local information
  mat.nbsa = matrix->n_nonzeros;
  // all this shit needs copying :(
  alloc1D_D(&(mat.sa), mat.nbsa);
  alloc1D_D(&(mat.jsa), mat.nbsa);
  for (int i = 0; i < mat.nbsa; ++i) {
    setElem1D_D(&(mat.sa), i+1, matrix->nnz[i]);
    setElem1D_I(&(mat.jsa), i+1, matrix->index_j[i]);
  }
  alloc1D_D(&(mat.isa), mat.n+1);
  for (int i = 0; i < mat.n+1; ++i) {
    setElem1D_I(&(mat.isa), i+1, matrix->edge_i[i]);
  }

  // create the RHS
  alloc2D_D(&f, mat.n, f_col);      // allocates the memory
  for (int row = 0; row < partition_height; ++row) {
    for (int col = 0; col < f_col; ++col) {
      // set the value
      setElem2D_D(&f, row, col, RHS[row * f_col + col]);
    }
  }
  
  //solution = &f;      // I guess .. but that doesn't work.

}

/// \brief Solves the system
void intel_SPIKE::solve_full() {

  spike(&pspike, &mat, &f, &info);
}
