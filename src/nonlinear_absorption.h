#include "nonlinear_response.h"

class NonlinearAbsorption : public NonlinearResponse {
public:
  NonlinearAbsorption(double Ui, double density_of_neutrals, double pressure, double fraction, double scaling, std::unique_ptr<Ionization::Ionization> rate)
    :Ui(Ui), density_of_neutrals(density_of_neutrals), pressure(pressure),
     fraction(fraction), scaling(scaling), rate(std::move(rate)) {}

  void calculate_temporal_response(const Radial& electric_field, const Array2D<double>& electron_density, Radial& response) override {
    for (int i = 0; i < response.Nradius; ++i) {
      for (int j = 0; j < response.Ntime; ++j) {
        const double E = electric_field.rt(i, j).real();
        const double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(E, 2);
        const double W = scaling * rate->ionization_rate(2*I);
        response.rt(i, j) = W / (I+1) * Ui * (fraction*density_of_neutrals*pressure - electron_density(i, j)) * Constants::epsilon_0 * Constants::c * E;
      }
    }
  }
  
  double Ui, density_of_neutrals, pressure, fraction, scaling;
  std::unique_ptr<Ionization::Ionization> rate;
};
