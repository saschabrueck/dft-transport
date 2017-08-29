/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_
#define CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_

#include "Types.H"
#include <vector>

namespace quadrature_types {
enum quadrature_type {
  NONE=0,
  CCGL=1,   // Complex Contour Gauss-Legendre
  GL=2,     // Gauss-Legendre
  GC=3,     // Gauss-Chebychev
  TS=4,     // Tanh-sinh
  TR=5,     // Trapezoidal Rule
  CCMR=6,   // Complex Contour Midpoint Rule
  MR=7,     // Midpoint Rule
};
} // namespace

/*! \class Quadrature
 *  
 *  \brief Provides access to abscissae and weights for quadratures
 *
 *  The class implements quadrature types by providing access to abscissae and
 *  weights for quadrature
 */
class Quadrature {
 public:
  std::vector<CPX> abscissae;
  std::vector<CPX> weights;
  Quadrature(quadrature_types::quadrature_type, double, double, int);

};

#endif
