#ifndef LINEAR_H_
#define LINEAR_H_

#include <functional>
#include <complex>

#include "medium.h"

class Linear {
public:
  Linear(std::function<double(double)> linear_index)
    :n(linear_index) {}

  virtual std::complex<double> kz(double kperp, double omega) const = 0;

  virtual double group_velocity(double omega0) const {
    double domega = 0.01 * omega0;
    double dn = (n(omega0 + domega) - n(omega0 - domega)) / (2 * domega);
    double k0 = (dn * omega0 + n(omega0)) / Constants::c;
    double vg = 1.0 / k0;
    return vg;
  }

protected:
  std::function<double(double)> n;
};


class FreeSpace : public Linear {
public:
  FreeSpace(std::function<double(double)> linear_index)
    :Linear(linear_index) {}
  
  std::complex<double> kz(double kperp, double omega) const override {
    double k = n(omega) * omega / Constants::c;
    return std::sqrt(std::pow(k, 2) - std::pow(kperp, 2));
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

  virtual double group_velocity(double omega0) const override {
    double np = Medium::pressurize(pressure, n, omega0);
    double domega = 0.01 * omega0;
    double np_plus = Medium::pressurize(pressure, n, omega0+domega);
    double np_minus = Medium::pressurize(pressure, n, omega0-domega);
    double dn = (np_plus - np_minus) / (2 * domega);

    double k0 = (dn * omega0 + np) / Constants::c;
    double vg = 1.0 / k0;
    return vg;
  }

private:
  double R, nclad, pressure;
};





#endif // LINEAR_H_
