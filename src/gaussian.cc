#include "gaussian.h"
#include "constants.h"
#include <iostream>

namespace Field {

  Gaussian::Gaussian(double wavelength, double waist, double focus, double length,
		     double phase, double delay, double energy,
                     double chirp, double gvd)
    :wavelength(wavelength), waist(waist), focus(focus),
     tau(length*1.699/2), phase(phase),
     delay(delay), energy(energy),
     chirp(chirp), gvd(gvd) {
    k0 = 2*Constants::pi / wavelength;
    omega0 = k0 * Constants::c;
    zr = k0 * std::pow(waist, 2) / 2;
    df = focus / (1 + std::pow(focus/zr, 2));
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

  std::complex<double> Gaussian::operator()(double distance, double r, double t) const {
    t -= delay;
    double z = distance - Constants::c * t;
    double w = radius(z);
    double C = curvature(z);
    double psi = gouy(z);
    std::complex<double> i(0, 1);
    std::complex<double> arg(0, 0);
    arg += -std::pow(r/w, 2);
    arg += 0.5*i*k0*std::pow(r, 2)*C;
    arg += -i*psi;
    
    double zds = std::pow(tau, 2) / (2*gvd);
    double T = tau * std::sqrt(std::pow(1+chirp*z/zds, 2) + std::pow(z/zds, 2));
    double eff_chirp = chirp + (1 + std::pow(chirp, 2)) * z/zds;
    double phi = std::atan(((1 + std::pow(chirp, 2)) * z + chirp) / zds);
    arg += -i*phi;
    
    if (focus == 0.0) {
      arg += -std::pow(t/T, 2) * (1.0 + i*(eff_chirp));
    }
    else {
      // double t_shifted = t + std::pow(r, 2) / (2*Constants::c*focus);
      double t_shifted = t + std::pow(r, 2) / (2*Constants::c*(focus-z));
      arg += -std::pow(t_shifted/T, 2) * (1.0 + i*eff_chirp);
    }

    arg += -i*omega0*t;
    arg += i*phase;

    double I0 = energy / (std::pow(w, 2) * std::pow(Constants::pi/2, 1.5) * T);
    double E0 = std::sqrt(2*I0 / (Constants::epsilon_0*Constants::c));
    return E0 * std::exp(arg);
  }

}
