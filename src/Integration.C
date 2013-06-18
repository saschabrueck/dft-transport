#include <exception>
#include <iostream>   // TODO: remove that later

#include "Integration.H"
#include "Integration_GL.H"

#ifndef PI
#define PI 3.141592653589793238462
#endif

/** \brief Initializer
 *
 *  Currently implemented methods for quadrature:
 *    - 'integral_type::CCGL'
 *      Complex contour Gauss Legendre. Abscissae and weights are
 *      precalculated up to machine precision, maximum number of
 *      abscissae is 120.
 *
 * \param integral_type The type of integral for which the abscissae/weights are to be loaded or calculated. See above for a list of supported values.
 * \param start Lower bound of range with nonzero state density
 * \param end Upper bound of range with nonzero state density
 * \param nabsc How many abscissae to use for the quadrature
 */
IntAbsc::IntAbsc(integral_type type, double start, double end, unsigned int nabsc)
{

  itype = type;

  if( start < end ) {
    band_start=start;
    band_end=end;
  }
  else if (start < end) {
    band_start=end;
    band_end=start;
  }
  //else {
  //  throw exc_intabsc;
  //}

  switch (itype) {
    case integral_type::CCGL:
      //if (nabsc >= GaussLegendre::absc.size() || nabsc == 0) {
      //  throw exc_intabsc;
      //}
      auto radius=(band_end-band_start)/2.0;
      auto center=(band_start+band_end)/2.0;
      const auto absc_prec=GaussLegendre::absc[nabsc-1];
      const auto weights_prec=GaussLegendre::weights[nabsc-1];
      for(int i=0; i<=nabsc-1; ++i) {
        // Shift the precalculated abscissa from [-1 1] to [0, pi] interval
        auto abs=absc_prec[i]*(PI-0.0)/2.0+(PI+0.0)/2.0;
        // Account for length of parametrization and the derivative along
        // the contour
        weights.push_back(weights_prec[i]*(PI-0)/2.0*radius*CPX(sin(abs),cos(abs)));
        // Calculate point on the contour
        absc.push_back(radius*CPX(-cos(abs),sin(abs))+center);
      }
      break;
    // No default as enum class is strongly typed
  }
};

class excIntAbsc: public std::exception {
  virtual const char* what() const throw() {
    return "Generic fault in IntAbsc()";
  }
} exc_intabsc;
