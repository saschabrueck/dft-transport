/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/**  \brief Compute density and transmission for a given energy
 *
 *   \author Sascha A. Brueck
 */
#ifndef SemiSelfConsistent_H
#define SemiSelfConsistent_H

#include "FEMGrid.H"
#include "Poisson.H"
#include <vector>

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;
extern PARAM *parameter;
extern VOLTAGE *voltage;

int semiselfconsistent(cp2k_csr_interop_type,cp2k_csr_interop_type,cp2k_csr_interop_type*,cp2k_csr_interop_type*,std::vector<double>,std::vector<contact_type>,std::vector<int>,std::vector<int>,double,transport_parameters);

#endif
