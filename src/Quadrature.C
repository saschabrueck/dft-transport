#include <exception>
#include <iostream>   // TODO: remove that later

#include "Quadrature.H"

#ifndef PI
#define PI 3.141592653589793238462
#endif

class excQuadrature: public std::exception {
 public:
  excQuadrature(std::string input) : 
    errmsg("Error calling Quadrature(): " + input){}
  ~excQuadrature() throw() {}
  const char* what() const throw() {return errmsg.c_str();}
 private:
  std::string errmsg;
};

#include "Quadrature_GL.H"

/** \brief Initializer
 *
 *  Currently implemented methods for quadrature:
 *    - 'quadrature_type::GL'
 *      Real line Gauss-Legendre. Abscissae and weights are precalculated
 *      up to machine precision, maximum number of abscissae is 140.
 *    - 'quadrature_type::CCGL'
 *      Complex contour Gauss Legendre. Same as above but abscissae are
 *      located on a half circle in the upper complex plane.
 *    - 'quadrature_type::ANPS'
 *      Areshkin-Nikolic Pole Summation. Quadrature by replacing the
 *      fermi function with a suitable alternative to sum up over
 *      the residues of the poles. For more details, see
 *      'Electron density and transport in top-gated graphene nanoribbon
 *      devices: First-principles Green function algorithms for systems 
 *      containing a large number of atoms', Areshkin&Nikolic, 
 *      Phys. Rev. B. 81, 2010
 *    - 'quadrature_type::GC'
 *      Real line Gauss-Chebychev. Abscissae and weights are calculated
 *      on the fly. This method allows for a pole on the upper end of the
 *      integration domain.
 *
 * \param quadrature_type The type of quadrature for which the 
 *                        abscissae/weights are to be loaded or 
 *                        calculated. See above for a list of 
 *                        supported values.
 * \param start Lower bound of range with nonzero state density in [eV]
 * \param end Upper bound of range with nonzero state density in [eV]
 * \param T Temperature in [K]
 * \param Ef Fermi level in [eV]
 * \param num_abscissae For GL, CCGL and GC: How many abscissae to use 
 *                      for the quadrature. For ANPS: exponent specifying
 *                      the precision, base e. Set to n to specify a
 *                      precision of exp(-6).
 */
Quadrature::Quadrature(quadrature_types::quadrature_type type, double start, 
                       double end, double T, double Ef, 
                       unsigned int num_abscissae) {
  double k=K_BOLTZMANN;
  my_type = type;
  if (start < end) {
    band_start = start;
    band_end = end;
  } else if (start > end) {
    band_start = end;
    band_end = start;
  } else {
    throw excQuadrature("Invalid range specified");
  }
  if (num_abscissae == 0) {
    throw excQuadrature("Invalid number of abscissae specified");
  }
  switch (my_type) {
    case quadrature_types::CCGL: {
      if (num_abscissae == 0) {
        throw excQuadrature("Invalid number of abscissae");
      } 
      double radius = (band_end - band_start) / 2.0;
      double center = (band_start + band_end) / 2.0;
      const std::vector<double> abscissae_precalc = 
                                    GaussLegendre::get_abscissae(num_abscissae-1);
      const std::vector<double> weights_precalc = 
                                    GaussLegendre::get_weights(num_abscissae-1);
      CPX fermi(1.0, 0.0);
      for (uint i = 0; i <= num_abscissae - 1; ++i) {
        // Shift the precalculated abscissa from [-1 1] to [0, pi] interval
        double abscissa = abscissae_precalc[i] * (PI - 0.0) / 2.0 +
                        (PI + 0.0) / 2.0;
        // Calculate point on the contour
        CPX contour_point = radius * CPX(-cos(abscissa), sin(abscissa)) +
                             center;
        abscissae.push_back(contour_point);
        // Fermi function at given point:
        if (T != 0) {
          fermi = 1.0 / (exp((contour_point - Ef) / (k * T)) + 1.0);
        }
        // Account for length of parametrization, the fermi level, and the 
        // derivative along the contour 
        weights.push_back(weights_precalc[i] * (PI - 0) / 2.0 * radius *
                          fermi * CPX(sin(abscissa), cos(abscissa)));
      }
      break;
    }
    case quadrature_types::GL: {
      if (num_abscissae == 0) {
        throw excQuadrature("Invalid number of abscissae");
      }
      const std::vector<double> abscissae_precalc = 
                                  GaussLegendre::get_abscissae(num_abscissae-1);
      const std::vector<double> weights_precalc = 
                                  GaussLegendre::get_weights(num_abscissae-1);
      double fermi = 1.0;
      for (uint i = 0; i <= num_abscissae - 1; ++i) {
        double abscissa = abscissae_precalc[i] * (band_end - band_start) / 2.0 + 
                       (band_end + band_start) / 2.0;
        abscissae.push_back(abscissa);
        if (T != 0) {
          fermi = 1.0 / (exp((abscissa - Ef) / (k * T)) + 1.0);
        } else if (abscissa > Ef) {
          fermi = 0.0;
        }
        weights.push_back(weights_precalc[i] * (band_end - band_start) / 2.0 *
                          fermi);
      }
      break;
    }
    case quadrature_types::ANPS: {
      // TODO: think about how to extend the interface to allow for
      //       specification of the precision (as exponent p of e)
      // TODO: check eqn. 8 in the paper and print a warning if it doesn't hold
      // TODO: this method doesn't work for T=0
      //auto test_begin = 10;
      //auto test_end = 30;
      int precision_exponent = 30;
      double precision = exp(-precision_exponent);
      // Loop for minimizing the number of poles for a given precision
      //for (auto vertical_poles = test_begin; vertical_poles <= test_end; 
      //     ++vertical_poles) {
        int vertical_poles = 20;
        // NOTE: I assume T_real = T but in fact T_real would be another
        // optimization parameter like vertical_poles.
        double T_real = T;
        double mu_real = band_start - (precision_exponent * k * T_real);
        double mu_imaginary = 2 * vertical_poles * PI * k * T;
        double T_imaginary = mu_imaginary / (precision_exponent * k);
        CPX mu_imaginary_cpx(0.0, mu_imaginary);
        CPX T_imaginary_cpx(0.0, T_imaginary);
        CPX z, residue;
        int n;

        n=0;
        while (true) {
          // Construct poles on 'vertical' line through the physical fermi level
          CPX z(Ef, (2 * n + 1) * PI * k * T);
          residue = -k * T * (1.0 / (1.0 + exp((z - mu_imaginary_cpx) /
                             (k * T_imaginary_cpx)))); // fermi_imaginary
          if (norm(residue) > precision) {
            abscissae.push_back(z);
            weights.push_back(residue);
            ++n;
          } else {
            break;
          }
        }

        n=0;
        // DEBUG
          std::cout << "precision: " << precision << std::endl;
          std::cout << "z\t\t\tresidue\t\tfermi_std\t\tfermi_real\t\tfermi-diff\t\ttotal\tnorm(total)\n";
        while (true) {
          // Construct poles on 'horizontal' line at height my_imaginary
          // from Ef downwards
          CPX z(Ef - (2 * n + 1) * PI * k * T_imaginary, mu_imaginary);
          residue = -k * T_imaginary_cpx *
                         ((1.0 / (1.0 + exp((z - Ef) / (k * T)))) - // fermi
                         (1.0 / (1.0 + exp((z - mu_real) / (k * T_real)))));
                                                               // fermi_real
          // DEBUG
            if (real(z) <= band_start) {
              std::cout << "overshoot\n";
            }
            std::cout << z << "\t" << -k * T_imaginary_cpx << "\t" <<
                (1.0 / (1.0 + exp((z - Ef) / (k * T)))) << "\t" <<
                (1.0 / (1.0 + exp((z - mu_real) / (k * T_real)))) << "\t" <<
                (1.0 / (1.0 + exp((z - Ef) / (k * T)))) -
                (1.0 / (1.0 + exp((z - mu_real) / (k * T_real)))) << "\t" <<
                residue << "\t" << norm(residue) << std::endl;
          if (norm(residue) > precision) {
            abscissae.push_back(z);
            weights.push_back(residue);
            ++n;
          } else {
            break;
          }
        }

        n=0;
        // DEBUG
          std::cout << "plus-side\n";
        while (true) {
          // Construct poles on 'horizontal' line at height my_imaginary
          // from Ef upwards
          CPX z(Ef + (2 * n + 1) * PI * k * T_imaginary, mu_imaginary);
          residue = -k * T_imaginary_cpx *
                         ((1.0 / (1.0 + exp((z - Ef) / (k * T)))) -
                         (1.0 / (1.0 + exp((z - mu_real) / (k * T_real)))));
          // DEBUG
            if (real(z) >= band_end) {
              std::cout << "overshoot\n";
            }
            std::cout << z << "\t" << -k * T_imaginary_cpx << "\t" <<
                (1.0 / (1.0 + exp((z - Ef) / (k * T)))) << "\t" <<
                (1.0 / (1.0 + exp((z - mu_real) / (k * T_real)))) << "\t" <<
                (1.0 / (1.0 + exp((z - Ef) / (k * T)))) -
                (1.0 / (1.0 + exp((z - mu_real) / (k * T_real)))) << "\t" <<
                residue << "\t" << norm(residue) << std::endl;
          if (norm(residue) > precision) {
            abscissae.push_back(z);
            weights.push_back(residue);
            ++n;
          } else {
            break;
          }
        }

        n=0;
        while (true) {
          // Construct poles on 'vertical' line through mu_real
          CPX z(mu_real, (2 * n + 1) * PI * k * T);
          residue = -k * T_real * (1.0 / (1.0 + exp((z - mu_imaginary_cpx) /
                             (k * T_imaginary_cpx)))); // fermi_imaginary
          if (norm(residue) > precision) {
            abscissae.push_back(z);
            weights.push_back(residue);
            ++n;
          } else {
            break;
          }
        }

    } // case ANPS
    case quadrature_types::GC: {
      if (num_abscissae <= 1) {
        throw excQuadrature("Invalid number of abscissae");
      }
      double fermi = 1.0;
      for (uint n = 1; n <= num_abscissae; ++n) {
        double abscissa = (cos(PI * (2 * n - 1.0) / (2.0 * num_abscissae)) *
                        (band_end - band_start) / 2.0 +
                        (band_end + band_start) / 2.0);
        abscissae.push_back(abscissa);
        if (T != 0) {
          fermi = 1.0 / (exp((abscissa - Ef) / (k * T)) + 1.0);
        } else if (abscissa > Ef) {
          fermi = 0.0;
        }
        weights.push_back((PI / num_abscissae) * 
                          (band_end - band_start) / 2.0 * fermi);
      }
      break;
    }
    default: {
      throw excQuadrature("Invalid integral type requested");
    }
  }
};
