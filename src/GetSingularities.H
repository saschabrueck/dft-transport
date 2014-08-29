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
int Execute(TCSR<double>*,TCSR<double>*,double*,double*,int*);
void write_bandstructure(int);
void get_propagating(vector< vector<double> >&,const vector<CPX>&,int);
~Singularities();

/// Vector containing the singularity points
std::vector<double> energies;
/// Lower integration limit
double energy_gs;
/// Valence band edge
double energy_vb;
/// Conduction band edge
double energy_cb;


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

double **energies_matrix;
double **derivatives_matrix;
double **curvatures_matrix;

double min_parabola(int,double*,double*,int,double&);
int determine_velocities(CPX*,CPX*,double,double*,double*,double*);
double determine_fermi(double,double,int);
int eigen(CPX*,CPX*,CPX*,CPX*,CPX,double*);
void follow_band(int);

};

#endif
