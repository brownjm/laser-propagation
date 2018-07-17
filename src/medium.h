#ifndef MEDIUM_H_
#define MEDIUM_H_

#include <functional>
#include "constants.h"


namespace Medium {

  double omega_to_microns(double omega);
  double index_vacuum(double);

  // index of air at 15 C, 1 atm (Ciddor 1996)
  double index_air(double omega);
  
  // index of argon at 15 C, 1 atm
  double index_argon(double omega);
  
  const std::function<double(double)> select_linear_index(const std::string& name);

  // calculate index as a function of pressure using Lorentz-Lorent equation
  double pressurize(double pressure, std::function<double(double)> index, double omega);

}

#endif // MEDIUM_H_
