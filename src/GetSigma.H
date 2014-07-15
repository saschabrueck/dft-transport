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
int Cutout(TCSR<CPX>*,int,CPX,transport_methods::transport_method,c_transport_type,MPI_Comm);
void Distribute(TCSR<CPX>*,MPI_Comm);
int GetSigma(MPI_Comm);

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
CPX *lambdaprol;
/// Propagating phase factors
CPX *lambdapror;

private:

int ndof;

int master_rank;

int do_delete_H;
int do_delete_sigma;
int do_delete_spsigdist;
int do_delete_inj;
int do_delete_spainjdist;

CPX energyp;//ONLY DEBUGGING
int complexenergypoint;
int compute_inj;

int n_cells;

int contact;

int bandwidth;
double colzerothr;
double eps_limit;
double eps_decay;

TCSR<CPX> *H0;
TCSR<CPX> *H1;
TCSR<CPX> *H1t;

};

#endif
