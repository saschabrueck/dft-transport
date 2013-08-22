/**  \brief Compute Boundary Self Energies and Injection
 *
 *   \author Sascha A. Brueck
 */
#ifndef __BOUNDARYSELFENERGY
#define __BOUNDARYSELFENERGY

#include "CSR.H"

class BoundarySelfEnergy {
public:
BoundarySelfEnergy();
~BoundarySelfEnergy();

int GetSigma(TCSR<CPX>*,CPX,int,transport_methods::transport_method,c_transport_type);

/// Self Energy
CPX *sigmal;
/// Self Energy
CPX *sigmar;
/// Injection
CPX *injl;
/// Injection
CPX *injr;
/// Number of Propagating States
int n_propagating;

private:

};

#endif
