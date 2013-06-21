#ifndef CPK_OMENDFT_TRANSPORT_PARDISO_H_
#define CPK_OMENDFT_TRANSPORT_PARDISO_H_

#include "Types.H"
#include <complex>
#include <omp.h>

namespace Pardiso {

void sparse_inv(TCSR<double> *matrix);
void sparse_inv(TCSR<CPX> *matrix);

}

#endif
