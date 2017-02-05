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
Singularities(transport_parameters,std::vector<contact_type>);
int Execute(cp2k_csr_interop_type,cp2k_csr_interop_type);
double determine_charge(double,int,int);
double determine_fermi(double,int);
void write_bandstructure(int,int,int,int);
std::vector< std::vector< std::vector<double> > > get_propagating(const std::vector<CPX>&);
int DensityFromBS(cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,std::vector<double>);
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
int dothederivs;
double eps_singularities;
double eps_mu;
int n_k;
double Temp;
double evfac;
int n_mu;

int iam;
int k_rank;
MPI_Comm bs_comm;
int size_bs_comm;
int rank_bs_comm;
MPI_Comm equal_bs_rank_comm;
std::vector<int> master_ranks;

std::vector<contact_type> contactvec;

std::vector< std::vector<double> > curvatures_extremum;
std::vector< std::vector<double> > kval_extremum;

std::vector< std::vector<double> > energies_matrix;
std::vector< std::vector<double> > derivatives_matrix;
std::vector< std::vector<double> > curvatures_matrix;

int determine_velocities(TCSR<double>**,TCSR<double>**,double,double*,double*,double*,int,int);
int determine_imaginary_bandstructure(TCSR<double>**,TCSR<double>**,double,CPX*,int);
void follow_band(int);
void add_full_to_scaled_cp2k_csr(cp2k_csr_interop_type&,double*,int,int,int,int,double,double);
int determine_density_from_bs(TCSR<double>**,TCSR<double>**,double**,double,double,int,int);

};

#endif
