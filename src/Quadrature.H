#ifndef CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_
#define CP2K_OMENDFT_TRANSPORT_QUADRATURE_H_

#include "Types.H"
#include <complex>

/** \enum quadrature_type
 *
 *  \brief User friendly access to quadrature types in \sa Quadrature
 *
 *  For a detailed description of the quadrature methods implemented, see
 *  \sa set_type
 *
 */
enum class quadrature_type
{
  CCGL
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
class Quadrature
{
  private:
    double band_start;
    double band_end;
    quadrature_type itype;
  public:
    /// Vector containing the abscissae
    std::vector<CPX> absc;
    
    /// Vector containing the corresponding weights
    std::vector<CPX> weights;

    /// Constructor for direct initialization
    Quadrature(quadrature_type, double, double, unsigned int);

    /// Get start of band
    double get_start() {return band_start;};

    /// Get end of band
    double get_end() {return band_end;};
};

/// \exception Generic fault in Quadrature
class exc_quadrature;

#endif
