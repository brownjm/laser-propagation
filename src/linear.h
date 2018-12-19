#ifndef LINEAR_H_
#define LINEAR_H_

#include <functional>
#include <complex>
#include <iostream>

#include "medium.h"

class Linear {
public:
  Linear(std::function<double(double)> linear_index)
    :n(linear_index) {}

  virtual std::complex<double> kz(double kperp, double omega) const = 0;

  virtual double group_velocity(double kperp, double omega) const {
    double domega = 1e-3 * omega;
    double kz_plus = kz(kperp, omega+domega).real();
    double kz_minus = kz(kperp, omega-domega).real();
    double vg = (2*domega) / (kz_plus - kz_minus);
    return vg;
  }


  std::function<double(double)> n;
};


class FreeSpace : public Linear {
public:
  FreeSpace(std::function<double(double)> linear_index)
    :Linear(linear_index) {}
  
  std::complex<double> kz(double kperp, double omega) const override {
    double k = n(omega) * omega / Constants::c;
    double k2 = std::pow(k, 2);
    double kperp2 = std::pow(kperp, 2);
    if (kperp2 > k2) {
      return 0.0;
    }
    else {
      double kzvalue = std::sqrt(std::pow(k, 2) - std::pow(kperp, 2));
      return kzvalue;
    }
  }
};



class Capillary : public Linear {
public:
  Capillary(std::function<double(double)> linear_index, double radius, double cladding_index,
	    double pressure)
    :Linear(linear_index), R(radius), nclad(cladding_index), pressure(pressure) {}
  
  std::complex<double> kz(double kperp, double omega) const override {
    double k0 = omega / Constants::c;
    double np = Medium::pressurize(pressure, n, omega);
    double k = np * k0;
    double beta = std::sqrt(std::pow(k, 2) - std::pow(kperp, 2));
    double eps = std::pow(nclad, 2); 
    double alpha = 1.0/(2*R) * std::pow(kperp / k0, 2) * (eps + 1) / std::sqrt(eps - 1);
    return std::complex<double>(beta, -alpha);
  }

private:
  double R, nclad, pressure;
};


#endif // LINEAR_H_
