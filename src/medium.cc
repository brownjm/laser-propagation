#include "medium.h"


double Medium::omega_to_microns(double omega) {
  return 1e6 * 2*Constants::pi*Constants::c / omega;
}
  
double Medium::index_vacuum(double) {
  return 1;
}

double Medium::index_air(double omega) {
  double microns = omega_to_microns(omega);
  double a = 0.05792105 / (238.0185 - std::pow(microns, -2));
  double b = 0.00167917 / (57.362 - std::pow(microns, -2));
  return 1 + a + b;
}

double Medium::index_argon(double omega) {
  double microns = omega_to_microns(omega);
  return 1 + 6.432135e-5 + 2.8606021e-2 / (144 - std::pow(microns, -2));
}
  
const std::function<double(double)> Medium::select_linear_index(const std::string& name) {
  if (name == "vacuum") return index_vacuum;
  else if (name == "argon") return index_argon;
  else if (name == "air") return index_air;
  else throw std::runtime_error("Unknown medium: " + name + "\n");
}

double Medium::pressurize(double pressure, std::function<double(double)> index, double omega) {
  double n0 = index(omega);
  double chi = (std::pow(n0, 2) - 1.0) / (std::pow(n0, 2) + 2.0);
  double np = std::sqrt((1.0 + 2.0*pressure*chi) / (1.0 - 2.0*chi));
  return np;
}




