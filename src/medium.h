#ifndef MEDIUM_H_
#define MEDIUM_H_

#include <complex>
#include <functional>
#include <gsl/gsl_spline.h>
#include "constants.h"
#include "io.h"

namespace Medium {

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
    Interpolated(const std::string& filename) {
      IO::read(filename, omega, index, absorption);
      omega_min = omega.front();
      omega_max = omega.back();
      spline_index = gsl_spline_alloc(gsl_interp_linear, omega.size());
      spline_absorption = gsl_spline_alloc(gsl_interp_linear, omega.size());

      gsl_spline_init(spline_index, omega.data(), index.data(), omega.size());
      gsl_spline_init(spline_absorption, omega.data(), absorption.data(), omega.size());
      
      acc_index = gsl_interp_accel_alloc();
      acc_absorption = gsl_interp_accel_alloc();
    }

    Interpolated(const Interpolated& other) {
      // copy the data that will be interpolated
      omega = other.omega;
      index = other.index;
      absorption = other.absorption;
      omega_min = omega.front();
      omega_max = omega.back();

      // create new spline interpolations
      spline_index = gsl_spline_alloc(gsl_interp_linear, omega.size());
      spline_absorption = gsl_spline_alloc(gsl_interp_linear, omega.size());

      gsl_spline_init(spline_index, omega.data(), index.data(), omega.size());
      gsl_spline_init(spline_absorption, omega.data(), absorption.data(), omega.size());
      
      acc_index = gsl_interp_accel_alloc();
      acc_absorption = gsl_interp_accel_alloc();
    }
    
    ~Interpolated() {
      gsl_spline_free(spline_index);
      gsl_spline_free(spline_absorption);
      gsl_interp_accel_free(acc_index);
      gsl_interp_accel_free(acc_absorption);
    }
    
    std::complex<double> operator()(double omega) {
      if (omega < omega_min) omega = omega_min;
      if (omega > omega_max) omega = omega_max;
      double n = gsl_spline_eval(spline_index, omega, acc_index);
      n = index_ethanol(omega).real();
      double k = gsl_spline_eval(spline_absorption, omega, acc_absorption);
      return std::complex<double>(n, k);
    }

  private:
    double omega_min, omega_max;
    std::vector<double> omega, index, absorption;
    gsl_spline* spline_index;
    gsl_spline* spline_absorption;
    gsl_interp_accel* acc_index;
    gsl_interp_accel* acc_absorption;
  };
  
  const std::function<std::complex<double>(double)> select_linear_index(const std::string& name);

  // calculate index as a function of pressure using Lorentz-Lorent equation
  std::complex<double> pressurize(double pressure, std::function<std::complex<double>(double)> index, double omega);

}

#endif // MEDIUM_H_
