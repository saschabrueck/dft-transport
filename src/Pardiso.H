#ifndef CPK_OMENDFT_TRANSPORT_PARDISO_H_
#define CPK_OMENDFT_TRANSPORT_PARDISO_H_

#include "Types.H"
#include <complex>
#include "CSR.H"
// will need to make some includes

namespace Pardiso {

void sparse_inv(TCSR<double> *matrix);
void sparse_inv(TCSR<CPX> *matrix);

}

#endif
