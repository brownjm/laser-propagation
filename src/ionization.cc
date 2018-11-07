#include "ionization.h"
#include "io.h"
#include "radial.h"
#include "constants.h"

namespace Ionization {

  Tabulated::Tabulated(const std::string& filename, double density_of_neutrals, double ionizing_fraction, double pressure, double scaling)
    :filename(filename), density_of_neutrals(density_of_neutrals),
     ionizing_fraction(ionizing_fraction), pressure(pressure), scaling(scaling) {
    std::vector<double> intensities, rates;
    IO::read(filename, intensities, rates);
    interp = std::make_unique<Interpolate>(intensities, rates);
  }

  void Tabulated::calculate_electron_density(const Radial& electric_field,
                                             Array2D<double>& electron_density) {
    const double dt = electric_field.time[1] - electric_field.time[0];
    for (int i = 0; i < electric_field.Nradius; ++i) {
      Util::IntegratorTrapz integrator(dt);
      for (int j = 0; j < electric_field.Ntime; ++j) {
        double E = electric_field.rt(i, j).real();
        double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(E, 2);
        double rate = scaling * interp->operator()(2*I);
        double probability = integrator.add(rate);
        electron_density(i, j) = density_of_neutrals * pressure * ionizing_fraction * probability;
      }
    }
  }
  
}
