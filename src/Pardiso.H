/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef CP2K_OMENDFT_TRANSPORT_PARDISO_H_
#define CP2K_OMENDFT_TRANSPORT_PARDISO_H_

#include <omp.h>
#include "Types.H"
#include "CSR.H"

// Support for various pardiso versions
#ifdef SRC_PARDISO
#define OLD_PARDISO_INTERFACE
#endif
#ifdef MKL_PARDISO
#define OLD_PARDISO_INTERFACE
#endif


extern "C" {

#ifndef SRC_PARDISO
// The sources based pardiso doesn't have a pardisoinit function
#ifdef OLD_PARDISO_INTERFACE
  void fortran_name(pardisoinit, PARDISOINIT)(void*, int*, int*);
#else
  void fortran_name(pardisoinit, PARDISOINIT)(void*, int*, int*, int*, 
                                              double*, int*);
#endif
#endif /* SRC_PARDISO */

#ifdef OLD_PARDISO_INTERFACE
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

#ifndef OLD_PARDISO_INTERFACE
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
#ifdef SRC_PARDISO
  iparm[3] = 80;
  iparm[20] = 1;
#endif

  int maxfct = 1;       // Maximum number of numerical factorizations.
  int mnum = 1;         // Which factorization to use.
  int msglvl = 0;       // Don't print statistical information in file.
  int error = 0;        // Initialize error flag to zero.

  for (int i = 0; i < 64; ++i) {
      pt[i] = 0;
  }

#ifndef SRC_PARDISO
  // Initialize pardiso
#ifdef OLD_PARDISO_INTERFACE
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
#endif /* SRC_PARDISO */



  // Factorization
#ifndef SRC_PARDISO
  int phase = 11; // Newer pardiso: reordering and symbolic factorization
#else
  int phase = 12; // SRC pardiso: analysis & numerical factorization
#endif
  //printf("entering phase %d\n", phase);
#ifdef OLD_PARDISO_INTERFACE
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

#ifndef SRC_PARDISO
  // Numerical factorization (not needed for SRC pardiso)
  phase = 22;
  //printf("entering phase %d\n", phase);
#ifdef OLD_PARDISO_INTERFACE
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
#endif /* SRC_PARDISO */
 
  // Back substitution and iterative refinement 
  phase = 33;
  iparm[7] = 2;       // Max number of refinement steps
  //printf("entering phase %d\n", phase);
#ifdef OLD_PARDISO_INTERFACE
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
#ifdef OLD_PARDISO_INTERFACE
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

#ifndef OLD_PARDISO_INTERFACE
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
#endif /* OLD_PARDISO_INTERFACE */

} /* namespace Pardiso */

#endif
