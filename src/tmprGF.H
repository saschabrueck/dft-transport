#ifndef CP2K_OMENDFT_TRANSPORT_TMPRGF_H_
#define CP2K_OMENDFT_TRANSPORT_TMPRGF_H_

#include "Types.H"
#include <complex>
#include <vector>
#include "CSR.H"
#include "rGF.H"

namespace tmprGF {

void sparse_invert(TCSR<CPX> *matrix, std::vector<int> Bmax);

}

#endif
