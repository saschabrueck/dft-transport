/**  \brief Get the points where the DOS is singular and determine the integration range
 *
 *   \author Sascha A. Brueck
 */
#ifndef __SINGULARITIES
#define __SINGULARITIES

#include "CSR.H"
#undef min
#undef max
#include <vector>
#include <valarray>
#include "Utilities.H"
#include <limits>

class Singularities {
public:
Singularities(c_transport_type,int);
int Execute(TCSR<double>*,TCSR<double>*,int*);
double determine_fermi(double,int);
void write_bandstructure(int);
std::vector< std::vector< std::vector<double> > > get_propagating(const std::vector<CPX>&);
~Singularities();

/// Vector containing the singularity points
std::vector< std::vector<double> > energies_extremum;
/// Lower integration limit
double energy_gs;
/// Valence band edge
double energy_vb;
/// Conduction band edge
double energy_cb;

std::vector<double> energies_vb;
std::vector<double> energies_cb;

private:
int iam, nprocs;
int ndof;
double eps_singularities;
int n_k;
int n_cells;
int bandwidth;
int noccunitcell;
double Temp;
int n_mu;

std::vector< std::vector<double> > curvatures_extremum;
std::vector< std::vector<double> > kval_extremum;

std::vector< std::vector<double> > energies_matrix;
std::vector< std::vector<double> > derivatives_matrix;
std::vector< std::vector<double> > curvatures_matrix;

int determine_velocities(CPX*,CPX*,double,double*,double*,double*);
int determine_imaginary_bandstructure(CPX*,CPX*,double,CPX*);
int eigen(CPX*,CPX*,CPX*,CPX*,CPX,double*);
void follow_band(int);

};

#endif
