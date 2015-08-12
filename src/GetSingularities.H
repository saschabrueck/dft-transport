/**  \brief Get the points where the DOS is singular and determine the integration range
 *
 *   \author Sascha A. Brueck
 */
#ifndef __SINGULARITIES
#define __SINGULARITIES

#include "CSR.H"
#include <vector>

class Singularities {
public:
Singularities(transport_parameters*,contact_type*);
int Execute(TCSR<double>*,TCSR<double>*);
double determine_fermi(double,int);
void write_bandstructure(int);
std::vector< std::vector< std::vector<double> > > get_propagating(const std::vector<CPX>&);
~Singularities();

/// Vector containing the singularity points
std::vector< std::vector<double> > energies_extremum;
/// Lower integration limit
double energy_gs;
/// Valence band edge
std::vector<double> energies_vb;
/// Conduction band edge
std::vector<double> energies_cb;

private:
double eps_singularities;
int n_k;
double Temp;
int n_mu;

int iam;
int k_rank;
MPI_Comm bs_comm;
int size_bs_comm;
int rank_bs_comm;
MPI_Comm equal_bs_rank_comm;
std::vector<int> master_ranks;

contact_type *contactvec;

std::vector< std::vector<double> > curvatures_extremum;
std::vector< std::vector<double> > kval_extremum;

std::vector< std::vector<double> > energies_matrix;
std::vector< std::vector<double> > derivatives_matrix;
std::vector< std::vector<double> > curvatures_matrix;

int determine_velocities(TCSR<double>**,TCSR<double>**,double,double*,double*,double*,int,int);
int determine_imaginary_bandstructure(TCSR<double>**,TCSR<double>**,double,CPX*,int);
void follow_band(int);

};

#endif
