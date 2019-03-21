#include <complex>
#include "medium.h"

double Medium::omega_to_microns(double omega) {
  return 1e6 * 2*Constants::pi*Constants::c / omega;
}
  
std::complex<double> Medium::index_vacuum(double) {
  return 1;
}

std::complex<double> Medium::index_air(double omega) {
  double microns = omega_to_microns(omega);
  double a = 0.05792105 / (238.0185 - std::pow(microns, -2));
  double b = 0.00167917 / (57.362 - std::pow(microns, -2));
  return 1 + a + b;
}

std::complex<double> Medium::index_argon(double omega) {
  double microns = omega_to_microns(omega);
  return 1 + 6.432135e-5 + 2.8606021e-2 / (144 - std::pow(microns, -2));
}

std::complex<double> Medium::index_ethanol(double omega) {
  double microns = omega_to_microns(omega);
  std::complex<double> gamma(0, 0.7); // for removing singularity at 3 microns
  std::complex<double> A = 0.0165*std::pow(microns, 2) / (std::pow(microns, 2) - 9.08 - gamma);
  double B = 0.8268*std::pow(microns, 2) / (std::pow(microns, 2) - 0.01039);
  std::complex<double> n = std::sqrt(1.0 + A + B);
  return n;
}


const std::function<std::complex<double>(double)> Medium::select_linear_index(const std::string& name) {
  if (name.find(".dat") != std::string::npos) {
    return Interpolated(name);
  }
  
  if (name == "vacuum") return index_vacuum;
  else if (name == "argon") return index_argon;
  else if (name == "air") return index_air;
  else if (name == "ethanol") return index_ethanol;
  else throw std::runtime_error("Unknown medium: " + name + "\n");
}

std::complex<double> Medium::pressurize(double pressure, std::function<std::complex<double>(double)> index, double omega) {
  std::complex<double> n0 = index(omega);
  std::complex<double> chi = (std::pow(n0, 2) - 1.0) / (std::pow(n0, 2) + 2.0);
  std::complex<double> np = std::sqrt((1.0 + 2.0*pressure*chi) / (1.0 - pressure*chi));
  return np;
}




