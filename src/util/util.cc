#include "util.h"
#include "constants.h"
#include "../core/radial.h"
#include <algorithm>
//#include <iostream>


double Util::energy(const Radial& field) {
  std::vector<double> aux(field.Nradius);
  double dt = field.time[1] - field.time[0];

  // integrate each row in time using the trapezoidal rule for uniform grids:
  // F = dx [1/2 f(x_0) + f(x_1) + f(x_1) + ... + 1/2 f(x_{N-1})]
  for (int i = 0; i < field.Nradius; ++i) {
    for (int j = 0; j < field.Ntime; ++j) {
      double I = 0.5 * Constants::epsilon_0 * Constants::c * std::norm(field.temporal(i, j));
      if ((j == 0) || (j == field.Ntime -1)) I *= 0.5;
      aux[i] += I;
    }
    aux[i] *= dt;
  }

  const auto& radius = field.radius;
  

  // double integral = aux[0] * std::pow(radius[0], 2) / 4;
  // for (int i = 1; i < field.Nradius; ++i) {
  //   integral += (radius[i]*aux[i] + radius[i-1]*aux[i-1]) / 2 * (radius[i] - radius[i-1]);
  // }
  // integral *= 2*Constants::pi;

  double integral = 0;
  for (int i = field.Nradius-1; i > 0; --i) {
      integral += (radius[i]*aux[i] + radius[i-1]*aux[i-1]) / 2 * (radius[i] - radius.at(i-1));
  }
  integral += aux[0] * std::pow(radius[0], 2) / 4;
  integral *= 2*Constants::pi;
  
  return integral;
}

double Util::energy(const std::vector<double>& radius, const std::vector<double>& time,
		    const std::vector<std::complex<double>>& E) {

  int Ntime = time.size();
  int Nradius = radius.size();
  std::vector<double> aux(Nradius);
  double dt = time[1] - time[0];

  // integrate each row in time using the trapezoidal rule for uniform grids:
  // F = dx [1/2 f(x_0) + f(x_1) + f(x_1) + ... + 1/2 f(x_{N-1})]
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Ntime; ++j) {
      double I = 0.5 * Constants::epsilon_0 * Constants::c * std::norm(E[i*Ntime + j]);
      if ((j == 0) || (j == Ntime -1)) I *= 0.5;
      aux[i] += I;
    }
    aux[i] *= dt;
  }


  
  // integral r f(r) dr from 0 to radius[0]
  double integral = aux[0] * std::pow(radius[0], 2) / 2;
  for (int i = 1; i < Nradius; ++i) {
    integral += (radius[i]*aux[i] + radius[i-1]*aux[i-1]) / 2 * (radius[i] - radius[i-1]);
  }
  integral *= 2*Constants::pi;
  
  return integral;
}

double Util::max_intensity(const std::vector<std::complex<double>>& E) {
  double maxI = 0;
  for (auto e : E) {
    double I = 0.5 * Constants::epsilon_0 * Constants::c * std::norm(e);
    if (I > maxI) {
      maxI = I;
    }
  }
  return maxI;
}

double Util::max_intensity(const Radial& field) {
  double maxI = 0;
  for (auto e : field.temporal.vec()) {
    double I = 0.5 * Constants::epsilon_0 * Constants::c * std::norm(e);
    if (I > maxI) {
      maxI = I;
    }
  }
  return maxI;
}

double Util::max_density(const Array2D<double>& rho) {
  return *std::max_element(std::begin(rho.values), std::end(rho.values));
}

Util::IntegratorSimps::IntegratorSimps(double dt)
  :dt(dt), f0(0), f1(0), f2(0), f3(0), F(0),
   a(3.0/8.0), b(19.0/24.0), c(-5.0/24.0), d(1.0/24.0) {}

double Util::IntegratorSimps::add(double new_value) {
  f0 = new_value * dt;
  F += a*f0 + b*f1 + c*f2 + d*f3;
  f3 = f2;
  f2 = f1;
  f1 = f0;
  return F;
}

Util::IntegratorTrapz::IntegratorTrapz(double dt)
  : dt(dt), F(0), old_value(0) {}

double Util::IntegratorTrapz::add(double new_value) {
  F += 0.5 * (new_value + old_value) * dt;
  old_value = new_value;
  return F;
}
