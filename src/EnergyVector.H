/**  \brief Distribute Matrices, construct Energyvector and start Density in parallel
 *
 *   \author Sascha A. Brueck
 */
#ifndef __ENERGYVECTOR
#define __ENERGYVECTOR

#include "CSR.H"

class Energyvector {
public:
Energyvector();
int Execute(TCSR<double>*,TCSR<double>*,std::vector<double>,std::vector<contact_type>,transport_parameters*);
~Energyvector();

private:
int determine_energyvector(std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&,std::vector< std::vector<int> >&,TCSR<double>*,TCSR<double>*,std::vector<double>&,std::vector<contact_type>,transport_parameters*);
int read_real_axis_energies(std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&);
int add_real_axis_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&,const std::vector< std::vector<double> > &,transport_parameters*);
int add_cmpx_cont_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&,transport_parameters*);
int iam, nprocs;

};

#endif
