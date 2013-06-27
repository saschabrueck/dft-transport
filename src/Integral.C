#include "CSR.H"
#include "Types.H"
#include "Density.H"
#include "Quadrature.H"

namespace Integral {

/** \brief Calculate the integral of the density matrix
 *
 *  Parallel implementation of integration of the density matrix.
 *
 * \param S_split Pointer to the local part of the overlap matrix
 * \param KohnSham_split Pointer to the local part of the Kohn-Sham matrix
 * \param P_split Pointer to the local part of the density matrix
 * \param parameters The parameters contained in a ParameterStruct struct
 */
void integrate(TCSR<double> *S_split, TCSR<double> *KohnSham_split, 
               TCSR<double> *P_split, ParameterStruct parameters) {
  int my_rank, num_processes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  TCSR<double> *KohnSham = new TCSR<double>(KohnSham_split, MPI_COMM_WORLD);
  TCSR<double> *S = new TCSR<double>(S_split, MPI_COMM_WORLD);
  TCSR<CPX> *P = new TCSR<CPX>(S->size, S->n_nonzeros, S->findx);
  P->init_variable(P->nnz,P->n_nonzeros);
  // TODO: decide what integral method to use on what energy range
  unsigned int num_abscissae = 100;
  Quadrature integral = Quadrature(quadrature_type::CCGL, parameters.band_start,
                                   parameters.band_end, parameters.T, 
                                   parameters.Ef, num_abscissae);
  int my_abscissae_start = ...;
  int my_abscissae_end = ...;
  int method = 2;       // 1 means GF, 2 means WF, for complex contour use 1 later
  for (auto i = my_abscissa_start; i < my_abscissa_end; ++i) {
    auto abscissa = integral.abscissae[i];
    auto weight = integral.weights[i];
    density(KohnSham, S, P, abscissa, weight, method, parameters);
  }
  P->reducescatterconvert(P, MPI_COMM_WORLD);
};

} // namespace
