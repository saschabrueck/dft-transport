/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
int Cutout(TCSR<CPX>*,contact_type,CPX,transport_methods::transport_method_type,MPI_Comm);
void Distribute(TCSR<CPX>*,MPI_Comm);
int GetSigma(MPI_Comm,transport_parameters);
int GetSigmaEig(MPI_Comm,transport_parameters);
int GetSigmaInv(MPI_Comm,transport_parameters);

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
int n_dec;
int eigval_degeneracy;

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
