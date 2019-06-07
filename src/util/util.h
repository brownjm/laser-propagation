#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <complex>

#include "../core/array2d.h"

class Radial;

namespace Util {

  double energy(const Radial& field);
  
  double energy(const std::vector<double>& radius,
		const std::vector<double>& time,
		const std::vector<std::complex<double>>& E);

  double max_intensity(const Radial& field);
  double max_intensity(const std::vector<std::complex<double>>& E);

  double max_density(const Array2D<double>& rho);


  class IntegratorSimps {
  public:
    IntegratorSimps(double dt);
    double add(double value);
    
  private:
    double dt, f0, f1, f2, f3, F;
    double a, b, c, d;
  };

  class IntegratorTrapz {
  public:
    IntegratorTrapz(double dt);
    double add(double value);

  private:
    double dt, F, old_value;
  };
}

#endif // UTIL_H_
