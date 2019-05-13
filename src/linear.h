#ifndef LINEAR_H_
#define LINEAR_H_

#include <functional>
#include <complex>
#include <iostream>

#include "medium.h"

namespace Linear {
  
class Base {
public:
  Base(Medium::IndexFunction linear_index);
  virtual std::complex<double> kz(double kperp, double omega) const = 0;
  virtual double group_velocity(double kperp, double omega) const;
  virtual double gvd(double kperp, double omega) const;

  Medium::IndexFunction n;
};


class FreeSpace : public Base {
public:
  FreeSpace(Medium::IndexFunction linear_index)
    :Base(linear_index) {}
  
  std::complex<double> kz(double kperp, double omega) const override;
};


class DiffractionLess : public Base {
public:
  DiffractionLess(Medium::IndexFunction linear_index)
    :Base(linear_index) {}

  std::complex<double> kz(double, double omega) const override;
};


class Capillary : public Base {
public:
  Capillary(Medium::IndexFunction linear_index, double radius, double cladding_index,
	    double pressure)
    :Base(linear_index), R(radius), nclad(cladding_index), pressure(pressure) {}

  std::complex<double> kz(double kperp, double omega) const override;

private:
  double R, nclad, pressure;
};

} // end namespace Linear
#endif // LINEAR_H_
