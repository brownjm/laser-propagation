#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include <complex>

#include "field.h"

namespace Field {

  // Gaussian spatial/temporal pulse focused with a parabolic mirror
  class Gaussian : public Field {
  public:
    // waist is spatial 1/e^2 radius, duration is full width half maximum (FWHM)
    // focus is positive for converging beams
    Gaussian(double wavelength, double waist, double focus, double duration, 
	     double phase, double delay, double energy, double chirp=0);
    std::complex<double> operator()(double radius, double time) const override;
    
  private:
    double wavelength, waist, focus, tau, phase, delay, energy, chirp, k0, omega0, zr, df;
    
    double radius(double z) const;
    double curvature(double z) const;
    double gouy(double z) const;
  };
}


#endif // GAUSSIAN_H_
