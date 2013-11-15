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

int Set_vars(CPX,MPI_Comm,MPI_Comm);
void Cutout(TCSR<CPX>*,int,MPI_Comm,c_transport_type);
void Distribute(TCSR<CPX>*,c_transport_type,MPI_Comm);
int GetSigma(transport_methods::transport_method,c_transport_type,MPI_Comm);

/// Self Energy
CPX *sigmal;
/// Self Energy
CPX *sigmar;
/// Distributed Self Energy
TCSR<CPX> *spsigmaldist;
/// Distributed Self Energy
TCSR<CPX> *spsigmardist;
/// Injection
CPX *injl;
/// Injection
CPX *injr;
/// Number of Propagating States
int n_propagating;

private:

int ndof;

int master_rank;

int do_delete_H;
int do_delete_sigma;
int do_delete_spsigdist;
int do_delete_inj;

int complexenergypoint;

TCSR<CPX> *H0;
TCSR<CPX> *H1;
TCSR<CPX> *H1t;

};

#endif
