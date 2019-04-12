#include "nonlinear_response.h"

class NonlinearAbsorption : public NonlinearResponse {
public:
  NonlinearAbsorption(double ionization_potential, double density_of_neutrals,
                      double pressure, double fraction,
                      Array2D<double>& ionization_rate) 
    :ionization_potential(ionization_potential),
     density_of_neutrals(density_of_neutrals*pressure),
     fraction(fraction), ionization_rate(ionization_rate) {}

  void calculate(const std::vector<double>& radius,
                 const std::vector<double>& time,
                 const Array2D<std::complex<double>>& electric_field,
                 const Array2D<double>& electron_density,
                 Array2D<std::complex<double>>& response) override {
    for (std::size_t i = 0; i < radius.size(); ++i) {
      for (std::size_t j = 0; j < time.size(); ++j) {
        double E = electric_field(i, j).real();
        double W = ionization_rate(i, j);
        double I = 0.5 * Constants::epsilon_0 * Constants::c * std::pow(E, 2);
        response(i, j) += W / (I+1) * ionization_potential * (fraction*density_of_neutrals - electron_density(i, j)) * Constants::epsilon_0 * Constants::c * E;
      }
    }
  }

private:
  double ionization_potential, density_of_neutrals, fraction;
  const Array2D<double>& ionization_rate;
};
