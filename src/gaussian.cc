#include "gaussian.h"
#include "constants.h"
#include <iostream>

namespace Field {

  Gaussian::Gaussian(double wavelength, double waist, double focus, double length,
		     double phase, double delay, double energy, double chirp)
    :wavelength(wavelength), waist(waist), focus(-focus),
     tau(length*1.699/2), phase(phase),
     delay(delay), energy(energy), chirp(chirp) {
    k0 = 2*Constants::pi / wavelength;
    omega0 = k0 * Constants::c;
    zr = k0 * std::pow(waist, 2) / 2;
    df = -focus / (1 + std::pow(focus/zr, 2));
  }

  double Gaussian::radius(double z) const {
    double b = std::pow(z/zr, 2);
    if (focus == 0.0) {
      return waist * std::sqrt(1 + b);
    }
    else {
      double a = std::pow(1 - z/focus, 2);
      return waist * std::sqrt(a + b);
    }
  }

  double Gaussian::curvature(double z) const {
    if (focus == 0.0) {
      return z / (std::pow(z, 2) + std::pow(zr, 2));
    }
    else {
      double a = z - df;
      return a / (std::pow(a, 2) + df * (focus - df));
    }
  }

  double Gaussian::gouy(double z) const {
    if (focus == 0.0) {
      return std::atan(z/zr);
    }
    else {
      double a = z - df;
      double b = std::sqrt(focus * df - std::pow(df, 2));
      return std::atan(a / b);
    }
  }
  
  std::complex<double> Gaussian::operator()(double r, double t) const {
    t -= delay;
    double z = Constants::c * t;
    double w = radius(z);
    double C = curvature(z);
    double psi = gouy(z);
    std::complex<double> i(0, 1);
    std::complex<double> arg(0, 0);
    arg += -std::pow(r/w, 2);
    arg += 0.5*i*k0*std::pow(r, 2)*C;
    arg += -i*psi;
    if (focus == 0.0) {
          arg += -std::pow(t/tau, 2);
          arg += -i*chirp*std::pow(t, 2);
    }
    else {
      double t_shifted = t - std::pow(r, 2) / (2*Constants::c*focus);
      arg += -std::pow(t_shifted/tau, 2);
      arg += -i*chirp*std::pow(t_shifted, 2);
    }

    arg += i*omega0*t;
    arg += i*phase;
    // double A = std::pow(2/Constants::pi, 0.75) / (w * std::sqrt(tau));
    // return A* waist/w * std::exp(arg);
    // std::cout << "A0 = " << A * waist/w << "\n";
    // std::cout << "w/w0=" << waist / w << "\n";
    double I0 = energy / (std::pow(waist, 2) * std::pow(Constants::pi/2, 1.5) * tau);
    double E0 = std::sqrt(2*I0 / (Constants::epsilon_0*Constants::c));
    //std::cout << "E0 = " << E0 << "\n";
    return E0 * std::exp(arg);
  }

}
