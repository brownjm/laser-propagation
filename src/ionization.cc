#include <iostream>
#include <gsl/gsl_sf_dawson.h>
#include "ionization.h"
#include "io.h"
#include "radial.h"
#include "constants.h"


Ionization::Ionization(const std::string& filename, double density_of_neutrals,
                       double pressure, double ionizing_fraction)
  :density_of_neutrals(density_of_neutrals*pressure), ionizing_fraction(ionizing_fraction) {
  IO::read(filename, intensity_values, rate_values);
    
  // interpolation
  spline = gsl_spline_alloc(gsl_interp_linear, intensity_values.size());
  gsl_spline_init(spline, intensity_values.data(), rate_values.data(), intensity_values.size());
  acc = gsl_interp_accel_alloc();
}

Ionization::~Ionization() {
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline);
}

void Ionization::calculate_electron_density(const Radial& electric_field,
                                            Array2D<double>& ionization_rate,
                                            Array2D<double>& electron_density) {

  double dt = electric_field.time[1] - electric_field.time[0];
  double eta = ionizing_fraction * density_of_neutrals * dt / 2;
  for (int i = 0; i < electric_field.Nradius; ++i) {
    double E = electric_field.rt(i, 0).real();
    double I = Constants::epsilon_0*Constants::c * std::pow(E, 2);
    ionization_rate(i, 0) = rate(I);
    electron_density(i, 0) = eta * ionization_rate(i, 0);
    for (int j = 1; j < electric_field.Ntime; ++j) {
      double E = electric_field.rt(i, j).real();
      double I = Constants::epsilon_0*Constants::c * std::pow(E, 2);
      ionization_rate(i, j) = rate(I);
      double exp = std::exp(-(ionization_rate(i, j) + ionization_rate(i, j-1)) / 2 * dt);
      electron_density(i, j) = exp * (electron_density(i, j-1) + eta * ionization_rate(i, j-1)) + eta * ionization_rate(i, j);
    }
  }
}

GenerateRate::GenerateRate(double Ip, double wavelength)
  :Z(1), n(1), l(0), m(0) {
  omega = 2*Constants::pi*Constants::c / wavelength * au_time;
  E0 = Ip / au_energy;
  F0 = std::pow(2*E0, 1.5);
  nstar = Z / std::sqrt(2*E0);
  lstar = nstar - 1;

}

double GenerateRate::adk(double intensity) const {
  double F = std::sqrt(2*intensity / (Constants::c*Constants::epsilon_0)) / au_field;
  double w = E0 * std::sqrt(3*F/(Constants::pi*F0)) * C(nstar, lstar) * f(l, m);
  w *= std::pow(2*F0/F, 2*nstar-std::abs(m)-1);
  w *= std::exp(-2.0*F0 / (3.0*F));
  return w / au_time;
}

double GenerateRate::mpi(double intensity) const {
  double F = std::sqrt(2*intensity / (Constants::c*Constants::epsilon_0)) / au_field;
  double nu0 = E0 / omega;
  double intpart;
  std::modf(nu0 + 1, &intpart);
  double K = intpart;
  double w = E0 * std::sqrt(2) / Constants::pi * std::pow(4, 2*nstar-std::abs(m)-1);
  w *= C(nstar, lstar) * f(l, m) * std::pow(nu0, 2*nstar+2*K-std::abs(m)-1.5);
  w *= std::exp(2*K-nu0) * dawson(std::sqrt(2*(K-nu0)));
  w *= std::pow(F/F0, 2*K);
  return w / au_time;
}

double GenerateRate::tunnel(double intensity) const {
  double F = std::sqrt(2*intensity / (Constants::c*Constants::epsilon_0)) / au_field;
  double w = std::sqrt(3.0*F / (Constants::pi*F0)) * std::pow(2.0*F0/F, 2*n);
  w *= C(n, l) * std::pow(F / (2*F0), m+1) * f(l, m) * E0;
  w *= std::exp(-2.0*F0 / (3.0*F));
  return w / au_time;
}

double GenerateRate::yudin(double intensity) const {
  double F = std::sqrt(2*intensity / (Constants::c*Constants::epsilon_0)) / au_field;
  double gamma = omega / F * std::sqrt(2*E0);
  double Gm = std::pow(1 + std::pow(gamma, 2), 0.5*std::abs(m) + 0.75);
  double w = E0 * std::sqrt(3*alpha(gamma) / (2.0*std::pow(gamma, 3)));
  w *= C(nstar, lstar) * f(l, m) * Gm;
  w *= std::pow(2*F0/F, 2*nstar-std::abs(m)-1);
  w *= A(E0, omega, gamma);
  w *= std::exp(-2.0*F0 / (3.0*F) * g(gamma));
  return w / au_time;
}

double GenerateRate::ilkov(double intensity) const {
  double F = std::sqrt(2*intensity / (Constants::c*Constants::epsilon_0)) / au_field;
  double gamma = omega / F * std::sqrt(2*E0);
  double w = E0 * std::sqrt(6/Constants::pi) * C(nstar, lstar) * f(l, m);
  w *= std::pow(2*F0/F, 2*nstar-std::abs(m)-1.5);
  w *= std::pow(1 + std::pow(gamma, 2), 0.5*std::abs(m) + 0.75);
  w *= A(E0, omega, gamma);
  w *= std::exp(-2.0*F0 / (3.0*F) * g(gamma));
  return w / au_time;
}

int GenerateRate::factorial(int n) const {
  return std::round(std::tgamma(n+1));
}

double GenerateRate::f(double l, double m) const {
  double a = (2*l + 1) * factorial(l + std::abs(m));
  double b = std::pow(2, std::abs(m)) * factorial(std::abs(m)) * factorial(l - std::abs(m));
  return a / b;
}

double GenerateRate::C(double n, double l) const {
  double a = std::pow(2, 2*n);
  double b = n * std::tgamma(n+l+1) * std::tgamma(n-l);
  return a / b;
}

double GenerateRate::g(double gamma) const {
  double a = 3.0/(2.0*gamma);
  double b = (1.0 + 1.0/(2.0*std::pow(gamma, 2))) * std::asinh(gamma);
  double c = std::sqrt(1 + std::pow(gamma, 2)) / (2.0*gamma);
  return a*(b - c);
}

double GenerateRate::beta(double gamma) const {
  return 2.0*gamma / std::sqrt(1 + std::pow(gamma, 2));
}

double GenerateRate::alpha(double gamma) const {
  double a = std::asinh(gamma);
  double b = gamma / std::sqrt(1 + std::pow(gamma, 2));
  return 2*(a - b);
}

double GenerateRate::dawson(double x) const {
  return gsl_sf_dawson(x);
}

double GenerateRate::A(double E0, double omega, double gamma) const {
  double nu = E0 / omega * (1 + 1.0/(2.0*std::pow(gamma, 2)));

  // terms in sum decrease in value as S increases
  // sum from smallest term to largest
  double result = 0;
  for (int S = 2000; S >= 0; --S) {
    double intpart;
    std::modf(E0/omega + 1, &intpart);
    double kappa = intpart + S;
    if (kappa >= nu) {
      double expt = std::exp(-alpha(gamma) * (kappa - nu));
      double arg = std::sqrt(beta(gamma) * (kappa - nu));
      double wm = dawson(arg);
      result += expt * wm;
    }
  }
  double prefactor = 4/std::sqrt(3*Constants::pi) * std::pow(gamma, 2) / (1 + std::pow(gamma, 2));
  return prefactor * result;
}
