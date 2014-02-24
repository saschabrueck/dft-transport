/**  \brief Get the points where the DOS is singular and determine the integration range
 *
 *   \author Sascha A. Brueck
 */
#ifndef __SINGULARITIES
#define __SINGULARITIES

#include "CSR.H"
#include <vector>
#include <valarray>
#include "Utilities.H"

class Singularities {
public:
Singularities(c_transport_type);
int Execute(TCSR<double>*,TCSR<double>*,int,double*,double*,int*);
void get_propagating(vector< vector<double> >&,const vector<CPX>&);
void delete_matrices();
~Singularities();

/// Vector containing the singularity points
std::vector<double> energies;
/// Lower integration limit
double energy_gs;

private:
int ndof;
double eps_singularities;
int n_k;
int n_cells;
int bandwidth;
int noccunitcell;

double *energies_matrix;
double *derivatives_matrix;
double *curvatures_matrix;

double min_parabola(int,double*,double*,int,double&);
int determine_velocities(CPX*,CPX*,double,double*,double*,double*);
void determine_fermi(double&,double);
int eigen(CPX*,CPX*,CPX*,CPX*,CPX,double*);
void write_bandstructure(int);
void follow_band();

};

#endif
