#include "ionization.h"
#include "io.h"
#include "radial.h"
#include "constants.h"

namespace Ionization {


  TabulatedRate::TabulatedRate(const std::string& filename, double scaling)
    :filename(filename), scaling(scaling) {
    std::vector<double> intensity, rate;
    IO::read(filename, intensity, rate);
    interp = std::make_unique<Interpolate>(intensity, rate);
  }

  double TabulatedRate::ionization_rate(double electric_field) {
    const double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(electric_field, 2);
    return scaling * interp->operator()(2*I);
  }

  IonizationModel::IonizationModel(double density_of_neutrals, double ionizing_fraction,
                                   double pressure, Rate& rate, int Nradius, int Ntime)
    :density_of_neutrals(pressure*density_of_neutrals), ionizing_fraction(ionizing_fraction),
     rate(rate), cached_rate(Nradius, Ntime) {}

  void IonizationModel::calculate_electron_density(const Radial& electric_field,
                                              Array2D<double>& electron_density) {
    const double dt = electric_field.time[1] - electric_field.time[0];
    for (int i = 0; i < electric_field.Nradius; ++i) {
      Util::IntegratorTrapz integrator(dt);
      for (int j = 0; j < electric_field.Ntime; ++j) {
        double E = electric_field.rt(i, j).real();
        double W = rate.ionization_rate(E);
        cached_rate(i, j) = W;
        double probability = integrator.add(W);
        electron_density(i, j) = density_of_neutrals * ionizing_fraction * probability;
      }
    }
  }
  
}
