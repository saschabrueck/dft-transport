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
int Execute(TCSR<double>*,TCSR<double>*,int,double*,int*,double*,double*,double*,int,int*,c_transport_type);
~Energyvector();

private:
int determine_energyvector(std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&,std::vector< std::vector<int> >&,TCSR<double>*,TCSR<double>*,double*,int*,c_transport_type,int);
void add_real_axis_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&,const std::vector< std::vector<double> > &,c_transport_type,int);
void add_cmpx_cont_energies(double,double,double,std::vector<CPX>&,std::vector<CPX>&,std::vector<transport_methods::transport_method>&,c_transport_type,int);
int iam, nprocs;

};

#endif
