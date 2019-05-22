#ifndef IONIZATION_H_
#define IONIZATION_H_

#include <vector>
#include <memory>
#include "interpolate.h"
#include "array2d.h"
#include "util.h"

class Radial;

class Ionization {
public:
  Ionization(const std::string& filename, double density_of_neutrals, double pressure,
             double ionizing_fraction);

  ~Ionization();

  double rate(double I) {
    return gsl_spline_eval(spline, I, acc);
  }
    
  void calculate_electron_density(const Radial& electric_field,
                                  Array2D<double>& ionization_rate,
                                  Array2D<double>& electron_density);

  double density_of_neutrals, ionizing_fraction;
  std::vector<double> intensity_values, rate_values;
  gsl_spline* spline;
  gsl_interp_accel* acc;
};


class GenerateRate {
public:
  GenerateRate(double Ip, double wavelength);

  double adk(double intensity) const;
  double mpi(double intensity) const;
  double tunnel(double intensity) const;
  double yudin(double intensity) const;
  double ilkov(double intensity) const;

private:
  const double au_energy = 4.3597441775e-18;
  const double au_time = 2.41888432650516e-17;
  const double au_field = 5.1422065211e11;
  int Z, n, l, m;
  double omega, E0, F0, nstar, lstar;
  
  int factorial(int n) const;
  double f(double l, double m) const;
  double C(double n, double l) const;
  double g(double gamma) const;
  double beta(double gamma) const;
  double alpha(double gamma) const;
  double dawson(double x) const;
  double A(double E0, double omega, double gamma) const;
};

#endif // IONIZATION_H_
