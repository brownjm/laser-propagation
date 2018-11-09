#include "nonlinear_response.h"

class NonlinearAbsorption : public NonlinearResponse {
public:
  NonlinearAbsorption(double ionization_potential, double density_of_neutrals,
                      double pressure, double fraction,
                      std::shared_ptr<Ionization::IonizationModel> ioniz_model)
    :ionization_potential(ionization_potential),
     density_of_neutrals(density_of_neutrals*pressure),
     fraction(fraction), model(ioniz_model) {}

  void calculate_temporal_response(const Radial& electric_field,
                                   const Array2D<double>& electron_density,
                                   Radial& response) override {
    for (int i = 0; i < response.Nradius; ++i) {
      for (int j = 0; j < response.Ntime; ++j) {
        const double E = electric_field.rt(i, j).real();
        const double W = model->cached_rate(i, j);
        const double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(E, 2);
        response.rt(i, j) = W / (I+1) * ionization_potential * (fraction*density_of_neutrals - electron_density(i, j)) * Constants::epsilon_0 * Constants::c * E;
      }
    }
  }
  
  double ionization_potential, density_of_neutrals, fraction;
  std::shared_ptr<Ionization::IonizationModel> model;
};
