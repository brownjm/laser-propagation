#include "nonlinear_response.h"

class Plasma : public NonlinearResponse {
public:
  Plasma(double collision, double pressure)
    : collision_time(collision / pressure) {}


  void calculate_response(const std::vector<double>& radius,
                          const std::vector<double>& time,
                          const Array2D<std::complex<double>>& electric_field,
                          const Array2D<double>& electron_density,
                          Array2D<std::complex<double>>& response) override {
    // integrate dJ/dt + J(t)/tau_c = e^2/m * rho(t) E(t)
    // using the exponential time differencing method:
    // S.M. Cox, P.C. Matthews, J. Comp. Phys. 176, 430 (2002)
    double dt = time[1] - time[0];
    double exp = std::exp(-dt / collision_time);
    double eta = dt * std::pow(Constants::e, 2) / (2 * Constants::m_e);
    for (std::size_t i = 0; i < radius.size(); ++i) {
      response(i, 0) += eta * electron_density(i, 0) * electric_field(i, 0).real();
      for (std::size_t j = 1; j < time.size(); ++j) {
        response(i, j) += exp * (response(i, j-1) + eta*electron_density(i, j-1)*electric_field(i, j-1).real()) + eta*electron_density(i, j)*electric_field(i, j).real();
      }
    }
  }

private:
  double collision_time;
};
