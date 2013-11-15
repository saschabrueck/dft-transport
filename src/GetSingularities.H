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
Singularities(TCSR<double> *,TCSR<double> *,int,double *,double *,int *,c_transport_type);
~Singularities();

/// Vector containing the singularity points
std::vector<double> energies;
/// Lower integration limit
double energy_gs;

private:
int ndof;
double min_parabola(int,double*,double*,int,double&);
int determine_velocities(CPX*,CPX*,double,double*,double*,c_transport_type);
};

#endif
