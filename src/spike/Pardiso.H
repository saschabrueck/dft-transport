#ifndef CP2K_OMENDFT_TRANSPORT_PARDISO_H_
#define CP2K_OMENDFT_TRANSPORT_PARDISO_H_

#include "Types.H"
#include <complex>
#include "CSR.H"

namespace Pardiso {

// TODO: This whole thing should be a templated class

void sparse_invert(TCSR<double> *matrix);
void sparse_invert(TCSR<CPX> *matrix);

void sparse_solve(TCSR<CPX> *matrix, CPX *RHS, int RHS_col, CPX *X);
void sparse_solve(TCSR<double> *matrix, double *RHS, int RHS_col, double *X);


//int get_pardiso_mtype();
}

#endif
