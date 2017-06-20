/**  \brief Distribute Matrices, construct Energyvector and start Density in parallel
 *
 *   \author Sascha A. Brueck
 */
#ifndef __ENERGYVECTOR
#define __ENERGYVECTOR

#include "libcp2k.h"

namespace distribution_methods {
    enum distribution_method_type {
        NO_DISTRIBUTION,
        SPLITSOLVE_DISTRIBUTION,
        FLOOR_DISTRIBUTION,
        CEILING_DISTRIBUTION,
        MASTER_DISTRIBUTION,
    };
}

class Energyvector {
public:
Energyvector();
int Execute(cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,cp2k_csr_interop_type*,std::vector<double>&,std::vector<contact_type>,std::vector<int>,std::vector<int>,double*,double*,transport_parameters);
~Energyvector();

private:
int distribute_and_execute(std::vector<CPX>,std::vector<CPX>,std::vector<CPX>,std::vector<CPX>,std::vector<CPX>,std::vector< std::vector<int> >,distribution_methods::distribution_method_type,int,cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,cp2k_csr_interop_type*,std::vector<double>&,std::vector<contact_type>,std::vector<int>,std::vector<int>,double*,double*,transport_parameters);
std::vector<int> get_tsizes(distribution_methods::distribution_method_type,int,std::vector<int>,std::vector<int>,int,int);
int write_transmission_current(std::vector<CPX>,std::vector<CPX>,std::vector<double>,std::vector<double>,transport_parameters);
int determine_energyvector(std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,std::vector< std::vector<int> >&,cp2k_csr_interop_type,cp2k_csr_interop_type,std::vector<double>&,std::vector<contact_type>,transport_parameters);
int assign_real_axis_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,const std::vector< std::vector<double> > &,int,transport_parameters);
int assign_cmpx_cont_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,double,int);
int iam, nprocs;

};

#endif
