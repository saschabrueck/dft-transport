/**  \brief Get the points where the DOS is singular and determine the integration range
 *
 *   \author Sascha A. Brueck
 */
#ifndef __SINGULARITIES
#define __SINGULARITIES

#include "CSR.H"

class Singularities {
public:
Singularities(TCSR<double> *,TCSR<double> *,ParameterStruct);
~Singularities();

/// Number of singularity points
int n_energies;
/// Array containing the singularity points
double* energies;
/// Lower integration limit
double energy_gs;
/// Upper integration limit
double energy_vbe;

private:
int determine_velocities(TCSR<double>*,TCSR<double>*,double,double&,double&,double*,int&,ParameterStruct);
};

#endif
