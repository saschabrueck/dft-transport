/**  \brief Compute density and transmission for a given energy
 *
 *   \author Sascha A. Brueck
 */
#ifndef __DENSITY_H
#define __DENSITY_H

#include "CSR.H"

int density(TCSR<double> *,TCSR<double> *,TCSR<double> *,TCSR<double> *,CPX,CPX,transport_methods::transport_method_type,std::vector<double>,std::vector<contact_type>,std::vector<result_type>&,std::vector<int>,std::vector<int>,transport_parameters,MPI_Comm);
 
#endif
