#ifndef __DiagScaLapack
#define __DiagScaLapack

#include "ScaLapack.H"
#include "array_tools.H"
#include "CSR.H"
#include "Types.H"

int diagscalapack(TCSR<double>*, TCSR<double>*, c_transport_type);

#endif

