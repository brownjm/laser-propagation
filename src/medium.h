#ifndef MEDIUM_H_
#define MEDIUM_H_

#include <complex>
#include <functional>
#include "constants.h"

namespace Medium {

  double omega_to_microns(double omega);
  std::complex<double> index_vacuum(double);

  // index of air at 15 C, 1 atm (Ciddor 1996)
  std::complex<double> index_air(double omega);
  
  // index of argon at 15 C, 1 atm
  std::complex<double> index_argon(double omega);

  // index of ethanol at 15 C
  std::complex<double> index_ethanol(double omega);
  
  const std::function<std::complex<double>(double)> select_linear_index(const std::string& name);

  // calculate index as a function of pressure using Lorentz-Lorent equation
  std::complex<double> pressurize(double pressure, std::function<std::complex<double>(double)> index, double omega);

}

#endif // MEDIUM_H_
