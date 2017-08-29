/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
double determine_free_charge(double,int);
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
