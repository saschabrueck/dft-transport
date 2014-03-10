#ifndef CP2K_OMENDFT_TRANSPORT_PARDISO_H_
#define CP2K_OMENDFT_TRANSPORT_PARDISO_H_

#include <omp.h>
#include "Types.H"
#include "CSR.H"

extern "C" {

#ifdef MKL_PARDISO
  void fortran_name(pardisoinit, PARDISOINIT)(void*, int*, int*);
#else
  void fortran_name(pardisoinit, PARDISOINIT)(void*, int*, int*, int*, 
                                              double*, int*);
#endif

#ifdef MKL_PARDISO
  int fortran_name(pardiso, PARDISO)(void*, int*, int*, int*, int*, int*,
                                     void*, int*, int*, int*, int*, int*, int*,
                                     void*, void*, int*);
#else
  int fortran_name(pardiso, PARDISO)(void*, int*, int*, int*, int*, int*, 
                                     void*, int*, int*, int*, int*, int*, int*,
                                     void*, void*, int*, double*);
#endif

}


namespace Pardiso {

template <typename T>
void sparse_solve(TCSR<T>* A, T* b, int b_cols, T* x);

#ifndef MKL_PARDISO
template <typename T>
void sparse_invert(TCSR<T>* A);
#endif

template <typename T>
int get_matrix_type(TCSR<T>* A);

} /* namespace Pardiso */


// Template definitions
namespace Pardiso {

/** \brief Function to solve a sparse linear system
 * 
 *  \param[in]        A
 *                    The coefficient matrix in CSR format.
 *
 *  \param[in]        b
 *                    The right hand side in column major format.
 *
 *  \param[in]        b_cols
 *                    The number of columns in the right hand side, i.e. the
 *                    number of right hand sides.
 *
 *  \param[in|out]    x
 *                    The array to contain the solution of the system (same
 *                    size as b).
 *
 *  \pre              The given matrix is a non-singular matrix (complex or 
 *                    real) in CSR format, b is an array containing the right
 *                    hand side to solve for.
 *
 *  \post             x contains the solution of the system.
 */
template <typename T>
void sparse_solve(TCSR<T>* A, T* b, int b_cols, T* x) {

  if (A->findx != 1) {
      A->change_findx(1);
  }

  int   numthreads = omp_get_max_threads();
  int   n = A->size;
  int*  ia = A->edge_i;
  int*  ja = A->index_j;
  T*    a = A->nnz;

  // Pardiso parameters
  int   mtype = get_matrix_type(A);
  void* pt[64];

  int   iparm[64];
  for (int i = 0; i < 64; ++i) {
      iparm[i] = 0;
  }
  iparm[0] = 1;         // No solver default
  iparm[1] = 2;         // Fill-in reordering from METIS
  iparm[2] = numthreads;
  iparm[7] = 2;         // Max numbers of iterative refinement steps
  iparm[9] = 13;        // Perturb the pivot elements with 1E-13
  iparm[10] = 1;        // Use nonsymmetric permutation and scaling MPS

  int maxfct = 1;       // Maximum number of numerical factorizations.
  int mnum = 1;         // Which factorization to use.
  int msglvl = 0;       // Don't print statistical information in file.
  int error = 0;        // Initialize error flag to zero.

  for (int i = 0; i < 64; ++i) {
      pt[i] = 0;
  }

  // Initialize pardiso
#ifdef MKL_PARDISO
  //printf("***** Using old interface *****\n");
  fortran_name(pardisoinit,PARDISOINIT)(pt, &mtype, iparm);
#else
  //printf("***** Using new interface *****\n");
  double dparm[64];
  int solver = 0;       // Use the sparse direct solver
  fortran_name(pardisoinit,PARDISOINIT)(pt, &mtype, &solver, iparm, dparm,
                                        &error);
  switch (error) {
    case 0:
      break;
    case -10:
      printf("No license file found\n");
      exit(1);
      break;
    case -11:
      printf("License expired\n");
      exit(1);
      break;
    case -12:
      printf("Wrong username or hostname");
      exit(1);
      break;
    default:
      printf("Other error");
      exit(1);
      break;
  }
#endif



  // Reordering and Symbolic Factorization.
  int phase = 11;
  //printf("entering phase %d\n", phase);
#ifdef MKL_PARDISO
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase,
                                &n, a, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                NULL, NULL, &error);
#else
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase,
                                &n, a, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                NULL, NULL, &error, dparm);
#endif
  if (error != 0) {
      printf("\nERROR during symbolic factorization: %d", error);
      exit(error);
  }

  // Numerical factorization
  phase = 22;
  //printf("entering phase %d\n", phase);
#ifdef MKL_PARDISO
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase,
                                &n, a, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                NULL, NULL, &error);
#else
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase,
                                &n, a, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                NULL, NULL, &error, dparm);
#endif
  if (error != 0) {
      printf("\nERROR during numerical factorization: %d", error);
      exit(error);
  }
 
  // Back substitution and iterative refinement 
  phase = 33;
  iparm[7] = 2;       // Max number of refinement steps
  //printf("entering phase %d\n", phase);
#ifdef MKL_PARDISO
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase,
                                &n, a, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                b, x, &error);
#else
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase,
                                &n, a, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                b, x, &error, dparm);
#endif
  if (error != 0){
      printf("\nERROR during solution: %d", error);
      exit(error);
  }

  // Release of memory
  phase = -1;
  //printf("entering phase %d\n", phase);
#ifdef MKL_PARDISO
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase, &n,
                                NULL, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                NULL, NULL, &error);
#else
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase, &n,
                                NULL, ia, ja, NULL, &b_cols, iparm, &msglvl,
                                NULL, NULL, &error, dparm);
#endif
  if (error != 0){
      printf("\nERROR during memory release: %d", error);
      exit(error);
  }
}

#ifndef MKL_PARDISO
/** \brief Function to selectively invert a CSR matrix using Pardiso
 * 
 *  The sparsity pattern in the input is preserved in the output, i.e. 
 *  the inversion is only calculated where the input matrix was non-zero.
 *
 *  \param[in]          A
 *                      The matrix to be inverted (in-place)
 *
 *  \pre                The given matrix is a non-singular matrix (complex
 *                      or real) in CSR format
 *
 *  \post               The input matrix has been replaced by its inverse
 *
 *  \todo               The comments in this piece of code should reflect the
 *                      current state of the pardiso documentation.
 */
template <typename T>
void sparse_invert(TCSR<T> *A) {

  if (A->findx != 1) {
      A->change_findx(1);
  }

  int   numthreads = omp_get_max_threads();
  int   n = A->size;
  int*  ia = A->edge_i;
  int*  ja = A->index_j;
  T*    a = A->nnz;

  // Input parameters for pardiso (see pardiso documentation)
  int   mtype = get_matrix_type(A);
  void* pt[64];

  int   iparm[64];
  for (int i = 0; i < 64; ++i) {
      iparm[i] = 0;
  }
  iparm[0] = 1;
  iparm[1] = 2;
  iparm[2] = numthreads;
  iparm[9] = 8;
  iparm[20] = 1;
  iparm[23] = 1;
  iparm[24] = 1;

  int maxfct = 1;
  int mnum = 1;
  int solver = 0;
  int msglvl = 0;
  int error = 0;

  int nrhs = 0;
  double dparam[64];

  fortran_name(pardisoinit,PARDISOINIT)(pt,&mtype,&solver,iparm,
      dparam,&error);

  switch (error) {
    case 0:
      break;
    case -10:
      printf("No license file found\n");
      exit(1);
      break;
    case -11:
      printf("License expired\n");
      exit(1);
      break;
    case -12:
      printf("Wrong username or hostname");
      exit(1);
      break;
    default:
      printf("Other error");
      exit(1);
      break;
  }

  // Reordering and Symbolic Factorization.
  int phase = 11;
  //printf("entering phase %d\n", phase);
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase, &n,
               a, ia, ja, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error,
               dparam);
  if (error != 0) {
      printf("\nERROR during symbolic factorization: %d", error);
      exit(error);
  }

  // Numerical factorization
  phase = 22;
  //printf("entering phase %d\n", phase);
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase, &n,
               a, ia, ja, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error,
               dparam);
  if (error != 0) {
      printf("\nERROR during numerical factorization: %d", error);
      exit(error);
  }

  // Sparsity preserving inversion
  phase = -22;
  //printf("entering phase %d\n", phase);
  iparm[35] = 0;
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase, &n,
               a, ia, ja, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error,
               dparam);
  if (error != 0) {
      printf("\nERROR during inversion: %d", error);
      exit(error);
  }

  // Release of memory
  phase = -1;
  //printf("entering phase %d\n", phase);
  fortran_name(pardiso,PARDISO)(pt, &maxfct, &mnum, &mtype, &phase, &n,
               a, ia, ja, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error,
               dparam);
  if (error != 0){
      printf("\nERROR during memory release: %d", error);
      exit(error);
  }
}
#endif

} /* namespace Pardiso */

#endif
