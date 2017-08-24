/**  \brief Compute density and transmission for a given energy
 *
 *   \author Sascha A. Brueck
 */
#ifndef SemiSelfConsistent_H
#define SemiSelfConsistent_H

#include "FEMGrid.H"
#include "Poisson.H"

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;
extern PARAM *parameter;
extern VOLTAGE *voltage;

int semiselfconsistent(cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,cp2k_csr_interop_type*,std::vector<double>,std::vector<contact_type>,std::vector<int>,std::vector<int>,double,transport_parameters);

#endif
