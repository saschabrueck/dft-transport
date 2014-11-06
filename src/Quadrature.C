#include <math.h>
#include "Quadrature.H"
#include "Blas.H"

struct QUADRATURE_Exception{
    QUADRATURE_Exception(const int line,const char* file) {std::cerr<<"Error in line "<<line<<" of file "<<file<<std::endl;}
};

/** \brief Constructor
 *
 *  Currently implemented methods for quadrature:
 *
 *    - 'quadrature_type::GL'
 *      Real line Gauss-Legendre. Abscissae and weights are calculated by
 *      constructing a symmetric companion matrix whose eigenvalues are the 
 *      abscissa. The weights are derived from the first eigenvector.
 *    - 'quadrature_type::CCGL'
 *      Complex contour Gauss Legendre. Same as above but abscissae are
 *      located on a half circle in the upper complex plane.
 *    - 'quadrature_type::GC'
 *      Real line Gauss-Chebychev. Abscissae and weights are calculated
 *      on the fly. This method allows for a pole on the upper end of the
 *      integration domain.
 *
 * \param[in]           type
 *                      The type of quadrature for which the abscissae/weights
 *                      are to be loaded or calculated. See above for a list of
 *                      supported values.
 *
 * \param[in]           band_start     
 *                      Lower bound of range with nonzero state density in
 *                      [eV].
 *
 * \param[in]           band_end
 *                      Upper bound of range with nonzero state density in
 *                      [eV].
 *
 * \param[in]           num_abscissae
 *                      How many abscissae to use for the quadrature.
 */
Quadrature::Quadrature(quadrature_types::quadrature_type type, double band_start, 
                       double band_end, int num_abscissae)
{
  if (num_abscissae <= 0 || band_end <= band_start) {
    type = quadrature_types::NONE;
  }
  switch (type) {
    case quadrature_types::NONE: {
      break;
    }
    case quadrature_types::CCGL: {
      double radius = (band_end - band_start) / 2.0;
      double center = (band_start + band_end) / 2.0;

      // Companion matrix (tridiagonal, symmetric)
      double* offdiagonal = new double[num_abscissae - 1];
      for (int i = 0; i < num_abscissae - 1; ++i) {
        offdiagonal[i] = (i + 1) / sqrt(4.0 * ((i + 1) * (i + 1)) - 1);
      }
      double* diagonal = new double[num_abscissae];
      for (int i = 0; i < num_abscissae; ++i) {
        diagonal[i] = 0.0;
      }

      // Solve
      int status = 0;
      double* workspace = new double[2 * num_abscissae - 1];
      double* eigenvectors = new double[num_abscissae * num_abscissae];
      c_dsteqr('I', num_abscissae, diagonal, offdiagonal, eigenvectors,
               num_abscissae, workspace, &status);
      if (status != 0) {
        throw QUADRATURE_Exception(__LINE__,__FILE__);
      }

      delete[] offdiagonal;
      delete[] workspace;


      // Eigenvalues are unshifted abscissae
      std::vector<double> abscissae_precalc(num_abscissae);
      for (int i = 0; i < num_abscissae; ++i) {
        abscissae_precalc[i] = diagonal[i];
      }
      delete[] diagonal;

      // 2 * square of first element of each eigenvector are unshifted weights
      std::vector<double> weights_precalc(num_abscissae);
      for (int i = 0; i < num_abscissae; ++i) {
        double eig1 = eigenvectors[i * num_abscissae]; // column major
        weights_precalc[i] = 2 * eig1 * eig1;
      }
      delete[] eigenvectors;

      for (int i = 0; i < num_abscissae; ++i) {
        // Shift the precalculated abscissa from [-1 1] to [0, pi] interval
        double abscissa = abscissae_precalc[i] * (M_PI - 0.0) / 2.0 +
                        (M_PI + 0.0) / 2.0;
        // Calculate point on the contour
        CPX contour_point = radius * CPX(-cos(abscissa), sin(abscissa)) +
                             center;
        abscissae.push_back(contour_point);
        weights.push_back(weights_precalc[i] * (M_PI - 0) / 2.0 * radius *
                          CPX(sin(abscissa), cos(abscissa)));
      }
      break;
    }
    case quadrature_types::GL: {
      if (num_abscissae <= 1) {
        throw QUADRATURE_Exception(__LINE__,__FILE__);
      }

      // Companion matrix (tridiagonal, symmetric)
      double* offdiagonal = new double[num_abscissae - 1];
      for (int i = 0; i < num_abscissae - 1; ++i) {
        offdiagonal[i] = (i + 1) / sqrt(4.0 * ((i + 1) * (i + 1)) - 1);
      }
      double* diagonal = new double[num_abscissae];
      for (int i = 0; i < num_abscissae; ++i) {
        diagonal[i] = 0.0;
      }

      // Solve
      int status = 0;
      double* workspace = new double[2 * num_abscissae - 1];
      double* eigenvectors = new double[num_abscissae * num_abscissae];
      c_dsteqr('I', num_abscissae, diagonal, offdiagonal, eigenvectors,
               num_abscissae, workspace, &status);
      if (status != 0) {
        throw QUADRATURE_Exception(__LINE__,__FILE__);
      }

      delete[] offdiagonal;
      delete[] workspace;


      // Eigenvalues are unshifted abscissae
      std::vector<double> abscissae_precalc(num_abscissae);
      for (int i = 0; i < num_abscissae; ++i) {
        abscissae_precalc[i] = diagonal[i];
      }
      delete[] diagonal;

      // 2 * square of first element of each eigenvector are unshifted weights
      std::vector<double> weights_precalc(num_abscissae);
      for (int i = 0; i < num_abscissae; ++i) {
        double eig1 = eigenvectors[i * num_abscissae]; // column major
        weights_precalc[i] = 2 * eig1 * eig1;
      }
      delete[] eigenvectors;


      
      for (int i = 0; i < num_abscissae; ++i) {

        double abscissa = (abscissae_precalc[i] * (band_end - band_start) / 2.0) + 
                          (band_end + band_start) / 2.0;
        abscissae.push_back(abscissa);

        weights.push_back(weights_precalc[i] * (band_end - band_start) / 2.0);

      }
      break;
    }
    case quadrature_types::GC: {
      for (int n = 1; n <= num_abscissae; ++n) {
        double abscissa = (cos(M_PI * (2 * n - 1.0) / (2.0 * num_abscissae)) *
                        (band_start - band_end) / 2.0 +
                        (band_start + band_end) / 2.0);
        abscissae.push_back(abscissa);
        weights.push_back(sqrt((abscissa - band_start) * (band_end - abscissa)) *
                          M_PI / num_abscissae);
      }
      break;
    }
    case quadrature_types::TS: {
      double step = 3.0 * 2.0 / (num_abscissae - 1);
      for (int n = 1; n <= num_abscissae; ++n) {
        double tau = step / 2.0 * (2 * n - num_abscissae - 1);
        double abscissa = tanh(M_PI / 2.0 * sinh(tau)) *
                        (band_end - band_start) / 2.0 +
                        (band_end + band_start) / 2.0;
        abscissae.push_back(abscissa);
        double weight = M_PI / 2.0 * step * cosh(tau) /
                        (cosh( M_PI / 2.0 * sinh(tau)) *
                        cosh( M_PI / 2.0 * sinh(tau))) *
                        (band_end - band_start) / 2.0;
        weights.push_back(weight);
      }
      break;
    }
    case quadrature_types::TR: {
      if (num_abscissae <= 1) {
        throw QUADRATURE_Exception(__LINE__,__FILE__);
      }
      double step = (band_end - band_start) / (num_abscissae - 1);
      for (int n = 1; n <= num_abscissae; ++n) {
        abscissae.push_back(band_start + (n - 1) * step);
        if (n == 1 || n == num_abscissae) { 
            weights.push_back(step / 2.0);
        } else {
            weights.push_back(step);
        }
      }
      break;
    }
    case quadrature_types::MR: {
      double step = (band_end - band_start) / num_abscissae;
      for (int n = 0; n < num_abscissae; ++n) {
        abscissae.push_back(band_start + (n + 0.5) * step);
        weights.push_back(step);
      }
      break;
    }
    case quadrature_types::CCMR: {
      double radius = (band_end - band_start) / 2.0;
      double center = (band_end + band_start) / 2.0;
      for (int n = 0; n < num_abscissae; ++n) {
        CPX theta = radius * exp(CPX(0.0, M_PI * (n + 0.5) / num_abscissae));
        abscissae.push_back(center - theta);
        weights.push_back(CPX(0.0, M_PI / num_abscissae) * theta);
      }
      break;
    }
    default: {
        throw QUADRATURE_Exception(__LINE__,__FILE__);
    }
  }
}
