#ifndef CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_
#define CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_

#include "cxx11_check.H"
#include "Types.H"
#include <vector>
#include <complex>

// k in [eV/K]
#ifndef K_BOLTZMANN
#define K_BOLTZMANN 8.617332478e-5;
#endif


/** \enum quadrature_type
 *
 *  \brief User friendly access to quadrature types in \sa Quadrature
 *
 */
namespace quadrature_types {
enum quadrature_type {
  CCGL=1,   // Complex Contour Gauss-Legendre
  GL=2,     // Gauss-Legendre
  ANPS=3,   // Areshkin-Nikolic Pole Summation
  GC=4,     // Gauss-Chebychev
  //TR=5    // Trapezoidal Rule
};
} // namespace

/** \class Quadrature
 *  
 *  \brief Provides access to abscissae and weights for quadratures
 *
 *  The class implements quadrature types by providing access to abscissae and
 *  weights for quadrature
 *
 *  \author Mauro Calderara (originally)
 *  \version $Revision$
 *  \date $Date$
 */
class Quadrature {
 public:
  std::vector<CPX> abscissae;    ///< Vector containing the abscissae
  std::vector<CPX> weights; ///< Vector containing the corresponding weights
  Quadrature(quadrature_types::quadrature_type, double, double, double, 
      double, unsigned int);
  double get_start() {return band_start;}; ///< Get start of band
  double get_end() {return band_end;}; ///< Get end of band

 private:
  double band_start;
  double band_end;
  double fermi_level;
  double temperature;
  quadrature_types::quadrature_type my_type;
};

/// \exception Generic fault in Quadrature
class excQuadrature;

#endif
