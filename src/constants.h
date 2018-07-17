#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_const_mksa.h>
#include <cmath>

namespace Constants {
  const double pi = M_PI;
  const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
  const double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
  const double e = GSL_CONST_MKSA_ELECTRON_CHARGE;
  const double m_e = GSL_CONST_MKSA_MASS_ELECTRON;
}

#endif // CONSTANTS_H
