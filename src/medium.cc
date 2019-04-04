#include <complex>
#include "medium.h"

namespace Medium {

  double omega_to_microns(double omega) {
    return 1e6 * 2*Constants::pi*Constants::c / omega;
  }
  
  std::complex<double> index_vacuum(double) {
    return 1;
  }

  std::complex<double> index_air(double omega) {
    double microns = omega_to_microns(omega);
    double a = 0.05792105 / (238.0185 - std::pow(microns, -2));
    double b = 0.00167917 / (57.362 - std::pow(microns, -2));
    return 1 + a + b;
  }

  std::complex<double> index_argon(double omega) {
    double microns = omega_to_microns(omega);
    return 1 + 6.432135e-5 + 2.8606021e-2 / (144 - std::pow(microns, -2));
  }

  std::complex<double> index_ethanol(double omega) {
    double microns = omega_to_microns(omega);
    std::complex<double> gamma(0, 0.7); // for removing singularity at 3 microns
    std::complex<double> A = 0.0165*std::pow(microns, 2) / (std::pow(microns, 2) - 9.08 - gamma);
    double B = 0.8268*std::pow(microns, 2) / (std::pow(microns, 2) - 0.01039);
    std::complex<double> C = 0.002*std::pow(microns, 2) / (std::pow(microns, 2) - 90.6 - gamma);
    std::complex<double> n = std::sqrt(1.0 + A + B + C);
    return n;
  }

  Interpolated::Interpolated(const std::string& filename) {
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

  Interpolated::Interpolated(const Interpolated& other) {
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

  Interpolated::~Interpolated() {
    gsl_spline_free(spline_index);
    gsl_spline_free(spline_absorption);
    gsl_interp_accel_free(acc_index);
    gsl_interp_accel_free(acc_absorption);
  }

  std::complex<double> Interpolated::operator()(double omega) {
    if (omega < omega_min) omega = omega_min;
    if (omega > omega_max) omega = omega_max;
    double n = gsl_spline_eval(spline_index, omega, acc_index);
    n = index_ethanol(omega).real();
    double k = gsl_spline_eval(spline_absorption, omega, acc_absorption);
    return std::complex<double>(n, k);
  }



  const IndexFunction select_linear_index(const std::string& name) {
    if (name.find(".dat") != std::string::npos) {
      return Interpolated(name);
    }
  
    if (name == "vacuum") return index_vacuum;
    else if (name == "argon") return index_argon;
    else if (name == "air") return index_air;
    else if (name == "ethanol") return index_ethanol;
    else throw std::runtime_error("Unknown medium: " + name + "\n");
  }

  std::complex<double> pressurize(double pressure, std::function<std::complex<double>(double)> index, double omega) {
    std::complex<double> n0 = index(omega);
    std::complex<double> chi = (std::pow(n0, 2) - 1.0) / (std::pow(n0, 2) + 2.0);
    std::complex<double> np = std::sqrt((1.0 + 2.0*pressure*chi) / (1.0 - pressure*chi));
    return np;
  }


} // end namespace Medium

