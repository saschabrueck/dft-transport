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

void Deallocate_Sigma();
void Deallocate_Gamma();
void Deallocate_Injection();
void Finalize();
int Set_master(MPI_Comm,MPI_Comm);
int Cutout(TCSR<CPX>*,contact_type,CPX,transport_methods::transport_method,MPI_Comm);
void Distribute(TCSR<CPX>*,MPI_Comm);
int GetSigma(MPI_Comm,int,transport_parameters*);
int GetSigmaEig(MPI_Comm,int,transport_parameters*);
int GetSigmaInv(MPI_Comm,int,transport_parameters*);

/// Self Energy
CPX *sigma;
/// Gamma
CPX *gamma;
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
double rcond;

private:
void create_M_matrix(CPX*,int,TCSR<CPX>**,int,CPX);

int master_rank;

CPX energy;

int compute_inj;
int compute_gamma;

contact_type contact;

TCSR<CPX> *H0;
TCSR<CPX> *H1;
TCSR<CPX> *H1t;

TCSR<CPX> **H;

};

#endif
