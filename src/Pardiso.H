#ifndef CP2K_OMENDFT_TRANSPORT_PARDISO_H_
#define CP2K_OMENDFT_TRANSPORT_PARDISO_H_

#include "Types.H"
#include <complex>
#include "CSR.H"
// will need to make some includes

namespace Pardiso {

void sparse_invert(TCSR<double> *matrix);
void sparse_invert(TCSR<CPX> *matrix);

}

#endif
