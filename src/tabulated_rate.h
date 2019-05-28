#ifndef TABULATED_RATE_H_
#define TABULATED_RATE_H_

#include "ionization.h"
#include "interpolate.h"

class TabulatedRate : public Ionization {
public:
  TabulatedRate(const std::string& filename, double density_of_neutrals, double pressure,
             double ionizing_fraction);

  ~TabulatedRate();

  double rate(double I) {
    return gsl_spline_eval(spline, I, acc);
  }
    
  void calculate_electron_density(const Radial& electric_field,
                                  Array2D<double>& ionization_rate,
                                  Array2D<double>& electron_density) override;

  double density_of_neutrals, ionizing_fraction;
  std::vector<double> intensity_values, rate_values;
  gsl_spline* spline;
  gsl_interp_accel* acc;
};


#endif // TABULATED_RATE_H_
