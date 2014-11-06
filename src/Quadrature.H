#ifndef CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_
#define CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_

#include "Types.H"
#include <vector>

/** \enum quadrature_type
 *
 *  \brief User friendly access to quadrature types in \sa Quadrature
 *
 */
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

/** \class Quadrature
 *  
 *  \brief Provides access to abscissae and weights for quadratures
 *
 *  The class implements quadrature types by providing access to abscissae and
 *  weights for quadrature
 *
 *  \version $Revision$
 *  \date $Date$
 */
class Quadrature {
 public:
  std::vector<CPX> abscissae;    ///< Vector containing the abscissae
  std::vector<CPX> weights; ///< Vector containing the corresponding weights
  Quadrature(quadrature_types::quadrature_type, double, double, int);

};

#endif
