#ifndef CP2K_OMENDFT_TRANSPORT_INTEGRATION_H_
#define CP2K_OMENDFT_TRANSPORT_INTEGRATION_H_ 

#include "Types.H"
#include <complex>

/** \enum integral_type
 *
 *  \brief User friendly access to integral types in \sa IntegralAbsc
 *
 *  For a detailed description of the integration methods implemented, see
 *  \sa set_type
 *
 */
enum class integral_type
{
  CCGL
};

/** \class IntegralAbsc
 *  
 *  \brief Provides access to abscissae and weights for integrals/quadratures
 *
 *  The class implements integral types by providing access to abscissae and
 *  weights for quadrature
 *
 *  \author Mauro Calderara (originally)
 *  \version $Revision$
 *  \date $Date$
 */
class IntAbsc
{
  private:
    double band_start;
    double band_end;
    integral_type itype;
  public:
    /// Vector containing the abscissae
    std::vector<CPX> absc;
    
    /// Vector containing the corresponding weights
    std::vector<CPX> weights;

    /// Constructor for direct initialization
    IntAbsc(integral_type, double, double, unsigned int);

    /// Get start of band
    double get_start() {return band_start;};

    /// Get end of band
    double get_end() {return band_end;};
};

/// \exception Generic fault in IntAbsc
class exc_intabsc;

#endif
