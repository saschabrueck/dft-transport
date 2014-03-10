#include "Pardiso.H"

namespace Pardiso {

/** \brief Function that returns the pardiso code for the matrix type
 *
 *  \param[in]        A
 *                    Matrix in TCSR format
 */
template <>
int get_matrix_type(TCSR<double>* A) {
  return 11;          // pardiso code for general real matrix
};
template <>
int get_matrix_type(TCSR<CPX>* A) {
  return 13;          // pardiso code for general complex matrix
};

} /* namespace Pardiso */
