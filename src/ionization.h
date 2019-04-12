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


#endif // IONIZATION_H_
