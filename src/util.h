#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <complex>

#include "array2d.h"

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
    IntegratorSimps(double dt)
      :dt(dt), f0(0), f1(0), f2(0), f3(0), F(0),
       a(3.0/8.0), b(19.0/24.0), c(-5.0/24.0), d(1.0/24.0) {}
    double add(double value) {
      f0 = value * dt;
      F += a*f0 + b*f1 + c*f2 + d*f3;
      f3 = f2;
      f2 = f1;
      f1 = f0;
      return F;
    }
    
  private:
    double dt, f0, f1, f2, f3, F;
    double a, b, c, d;
  };

  class IntegratorTrapz {
  public:
    IntegratorTrapz(double dt)
      : dt(dt), F(0), old_value(0) {}

    double add(double new_value) {
      F += 0.5 * (new_value + old_value) * dt;
      old_value = new_value;
      return F;
    }

  private:
    double dt, F, old_value;
  };
}

#endif // UTIL_H_
