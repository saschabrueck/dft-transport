/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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

/*!  \brief Distribute matrices from CP2K to the energy points, call bandstructure calculator, construct energyvector and call density evaluation in parallel
 *
 *   \author Sascha A. Brueck
 */
class Energyvector {
public:
Energyvector();
int Execute(cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,cp2k_csr_interop_type*,std::vector<double>&,std::vector<contact_type>,std::vector<int>,std::vector<int>,double*,transport_parameters);
~Energyvector();

private:
int distribute_and_execute(std::vector<CPX>,std::vector<CPX>,std::vector<CPX>,std::vector<CPX>,std::vector<CPX>,std::vector< std::vector<int> >,distribution_methods::distribution_method_type,int,cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,cp2k_csr_interop_type*,std::vector<double>&,std::vector<contact_type>,std::vector<int>,std::vector<int>,double*,transport_parameters);
std::vector<int> get_tsizes(distribution_methods::distribution_method_type,int,std::vector<int>,std::vector<int>,int,int);
int write_transmission_current(std::vector<CPX>,std::vector<CPX>,std::vector<double>,std::vector<double>,transport_parameters);
int determine_energyvector(std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,std::vector< std::vector<int> >&,cp2k_csr_interop_type,cp2k_csr_interop_type,std::vector<double>&,std::vector<contact_type>,transport_parameters);
int assign_real_axis_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,const std::vector< std::vector<double> > &,int,transport_parameters);
int assign_cmpx_cont_energies(double,double,std::vector<CPX>&,std::vector<CPX>&,std::vector<CPX>&,double,int);
int iam, nprocs;

};

#endif
