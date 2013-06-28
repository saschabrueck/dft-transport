#include "Pardiso.H"
#include <cstddef>
#include <omp.h>


extern "C" {

    void fortran_name(pardisoinit,PARDISOINIT)(void*,int*,int*,int*,\
        double*,int*);
    int fortran_name(pardiso,PARDISO)(void*,int*,int*, int*,int*,int*,\
        CPX*,int*,int*,int*,int*,int*,int*,CPX*,CPX*,int*,double*);
}

namespace Pardiso {

// Use overload instead of template, I don't want to rely on RTTI

/** \brief Function to invert a CSR matrix using Pardiso
 * 
 *  This function inverts sparse matrices using the Pardiso library.
 *  The sparsity pattern in the input is preserved in the output, that
 *  the inversion is only calculated where the input matrix was non-zero.
 *
 *  \param *matrix The matrix to be inverted (in-place)
 *  \pre The given matrix is a non-singular matrix (complex or real) 
 *       in CSR format
 *  \post The input matrix has been replaced by its inverse
 */
void sparse_inv(TCSR<CPX> *matrix) {
  // Input parameters for pardiso (see pardiso documentation)
  int maxfct = 1;
  int mnum = 1;
  int mtype = 13;     // general linear complex
  int solver = 0;
  int n = matrix->size;
  int msglvl = 0;
  int error = 0;
  int phase = 0;
  int numthreads = omp_get_max_threads();
  int nrhs = 0;

  void* handle[64];
  int iparam[64];
  double dparam[64];

  fortran_name(pardisoinit,PARDISOINIT)(handle,&mtype,&solver,iparam,\
      dparam,&error);

  iparam[0] = 1;
  iparam[1] = 2;
  iparam[2] = numthreads;
  iparam[3] = 0;
  iparam[4] = 0;
  iparam[5] = 0;
  iparam[7] = 0;
  iparam[9] = 8;
  iparam[10] = 0;
  iparam[11] = 0;
  iparam[12] = 0;
  iparam[17] = 0;
  iparam[18] = 0;
  iparam[20] = 1;
  iparam[23] = 1;
  iparam[24] = 1;

  phase = 11;   // reorder and symbolic factorization
  fortran_name(pardiso,PARDISO)(handle,&maxfct,&mnum,&mtype,&phase,&n,
      matrix->nnz,matrix->edge_i,matrix->index_j,NULL,&nrhs,
      iparam,&msglvl,NULL,NULL,&error,dparam);

  phase = 22;   // factorize
  fortran_name(pardiso,PARDISO)(handle,&maxfct,&mnum,&mtype,&phase,&n,
      matrix->nnz,matrix->edge_i,matrix->index_j,NULL,&nrhs,
      iparam,&msglvl,NULL,NULL,&error,dparam);

  phase = -22;  // sparsity preserving inversion
  iparam[35] = 0;
  fortran_name(pardiso,PARDISO)(handle,&maxfct,&mnum,&mtype,&phase,&n,
      matrix->nnz,matrix->edge_i,matrix->index_j,NULL,&nrhs,
      iparam,&msglvl,NULL,NULL,&error,dparam);

}

} // namespace
