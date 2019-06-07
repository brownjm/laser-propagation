#ifndef MEDIUM_H_
#define MEDIUM_H_

#include <complex>
#include <functional>
#include <gsl/gsl_spline.h>
#include "../util/constants.h"
#include "../util/io.h"

namespace Medium {

  using IndexFunction = std::function<std::complex<double>(double)>;
    
  double omega_to_microns(double omega);
  std::complex<double> index_vacuum(double);

  // index of air at 15 C, 1 atm (Ciddor 1996)
  std::complex<double> index_air(double omega);
  
  // index of argon at 15 C, 1 atm
  std::complex<double> index_argon(double omega);

  // index of ethanol at 15 C
  std::complex<double> index_ethanol(double omega);

  // tabulated data from file
  class Interpolated {
  public:
    Interpolated(const std::string& filename);
    Interpolated(const Interpolated& other);
    ~Interpolated();
    std::complex<double> operator()(double omega);

  private:
    double omega_min, omega_max;
    std::vector<double> omega, index, absorption;
    gsl_spline* spline_index;
    gsl_spline* spline_absorption;
    gsl_interp_accel* acc_index;
    gsl_interp_accel* acc_absorption;
  };
  
  const IndexFunction select_linear_index(const std::string& name);

  // calculate index as a function of pressure using Lorentz-Lorent equation
  std::complex<double> pressurize(double pressure, IndexFunction index, double omega);

}

#endif // MEDIUM_H_
