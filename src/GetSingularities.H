/**  \brief Get the points where the DOS is singular and determine the integration range
 *
 *   \author Sascha A. Brueck
 */
#ifndef __SINGULARITIES
#define __SINGULARITIES

#include "CSR.H"

class Singularities {
public:
Singularities(TCSR<double> *,TCSR<double> *,c_transport_type);
~Singularities();

/// Vector containing the singularity points
std::vector<double> energies;
/// Lower integration limit
double energy_gs;
/// Upper integration limit
double energy_vbe;

private:
double min_parabola(int,double*,double*,int,double&);
int determine_velocities(TCSR<double>*,TCSR<double>*,double,int,double*,double*,c_transport_type);
};

#endif
