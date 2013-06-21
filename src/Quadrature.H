#ifndef CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_
#define CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_

#include "Types.H"
#include <complex>

/** \enum quadrature_type
 *
 *  \brief User friendly access to quadrature types in \sa Quadrature
 *
 *
 */
enum class quadrature_type {
  CCGL,   // Complex Contour Gauss Legendre
  GL,     // Gauss Legendre
  ANPS,   // Areshkin-Nikolic Pole Summation
  GC,     // Gauss Chebychev
  //TR    // Trapezoidal Rule
};

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
  Quadrature(quadrature_type, double, double, double, double, unsigned int);
  double get_start() {return band_start;}; ///< Get start of band
  double get_end() {return band_end;}; ///< Get end of band

 private:
  const double k = 8.617332478e-5;        // k in [eV/K]
  double band_start;
  double band_end;
  double fermi_level;
  double temperature;
  quadrature_type my_type;
};

/// \exception Generic fault in Quadrature
class excQuadrature;

#endif
