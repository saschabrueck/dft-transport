#ifndef SemiSelfConsistent_H
#define SemiSelfConsistent_H

#include "FEMGrid.H"
#include "Poisson.H"
#include "GetSingularities.H"
#include "EnergyVector.H"
#include <numeric>
#include <iterator>

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;
extern PARAM *parameter;
extern VOLTAGE *voltage;

int semiselfconsistent(TCSR<double>* , TCSR<double>* , transport_parameters* );

#endif
