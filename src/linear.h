#ifndef LINEAR_H_
#define LINEAR_H_

#include <functional>
#include <complex>
#include <iostream>

#include "medium.h"

class Linear {
public:
  Linear(std::function<std::complex<double>(double)> linear_index)
    :n(linear_index) {}

  virtual std::complex<double> kz(double kperp, double omega) const = 0;

  virtual double group_velocity(double kperp, double omega) const {
    double domega = 1e-3 * omega;
    double kz_plus = kz(kperp, omega+domega).real();
    double kz_minus = kz(kperp, omega-domega).real();
    double vg = (2*domega) / (kz_plus - kz_minus);
    return vg;
  }


  std::function<std::complex<double>(double)> n;
};


class FreeSpace : public Linear {
public:
  FreeSpace(std::function<std::complex<double>(double)> linear_index)
    :Linear(linear_index) {}
  
  std::complex<double> kz(double kperp, double omega) const override {
    std::complex<double> k = n(omega) * omega / Constants::c;
    std::complex<double> k2 = std::pow(k, 2);
    double kperp2 = std::pow(kperp, 2);
    if (kperp2 > std::real(k2)) {
      return 0.0;
    }
    else {
      std::complex<double> kzvalue = std::sqrt(std::pow(k, 2) - std::pow(kperp, 2));
      return kzvalue;
    }
  }
};



class Capillary : public Linear {
public:
  Capillary(std::function<std::complex<double>(double)> linear_index, double radius, double cladding_index,
	    double pressure)
    :Linear(linear_index), R(radius), nclad(cladding_index), pressure(pressure) {}
  
  std::complex<double> kz(double kperp, double omega) const override {
    double k0 = omega / Constants::c;
    std::complex<double> np = Medium::pressurize(pressure, n, omega);
    std::complex<double> k = np * k0;
    std::complex<double> beta = std::sqrt(std::pow(k, 2) - std::pow(kperp, 2));
    double eps = std::pow(nclad, 2); 
    double alpha = 1.0/(2*R) * std::pow(kperp / k0, 2) * (eps + 1) / std::sqrt(eps - 1);
    return beta + std::complex<double>(0, -alpha);
  }

private:
  double R, nclad, pressure;
};


#endif // LINEAR_H_
