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

int Set_master(MPI_Comm,MPI_Comm);
int Cutout(TCSR<CPX>*,contact_type,CPX,transport_methods::transport_method,MPI_Comm);
void Distribute(TCSR<CPX>*,MPI_Comm);
int GetSigma(MPI_Comm,int,transport_parameters*);

/// Self Energy
CPX *sigma;
/// Distributed Self Energy
TCSR<CPX> *spsigmadist;
/// Injection
CPX *inj;
/// Injection
TCSR<CPX> *spainjdist;
/// Number of Propagating States
int n_propagating;
/// Propagating phase factors
CPX *lambdapro;

private:
int master_rank;

CPX energy;

int compute_inj;

contact_type contact;

TCSR<CPX> *H0;
TCSR<CPX> *H1;
TCSR<CPX> *H1t;

TCSR<CPX> *spA;
TCSR<CPX> *spB;

};

#endif
