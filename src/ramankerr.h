#include "nonlinear_response.h"

class RamanKerr : public NonlinearResponse {
public:
  RamanKerr(double n2, double n0, double fraction, double gamma, double lambda, double pressure)
    :fraction(fraction), gamma(gamma), lambda(lambda) {
    chi3 = 4.0/3.0 * Constants::epsilon_0 * Constants::c * n2 * pressure * std::pow(n0, 2);
    R0 = (std::pow(gamma, 2)/4.0 + std::pow(lambda, 2)) / lambda;
  }

  void calculate(const std::vector<double>& radius,
                 const std::vector<double>& time,
                 const Array2D<std::complex<double>>& electric_field,
                 const Array2D<double>&,
                 Array2D<std::complex<double>>& response) override {
    double dt = time[1] - time[0];
    std::complex<double> imagi(0, 1);
    std::complex<double> exp = std::exp((-gamma/2.0 + imagi*lambda) * dt);
    double eta = R0 * dt / 2.0;
    for (std::size_t i = 0; i < radius.size(); ++i) {
      double E = electric_field(i, 0).real();
      double R = eta * std::pow(E, 2);
      response(i, 0) += Constants::epsilon_0 * chi3 * ((1 - fraction) * std::pow(E, 2)
                                                       + fraction * R) * E;
      double Rold = R;
      double Eold = E;
      for (std::size_t j = 1; j < time.size(); ++j) {
        E = electric_field(i, j).real();
        R = std::imag(exp*Rold + eta * (std::pow(E, 2) + eta * std::pow(Eold, 2)));
        response(i, j) += Constants::epsilon_0 * chi3 * ((1 - fraction) * std::pow(E, 2)
                                                         + fraction * R) * E;
        Eold = E;
        Rold = R;
      }
    }
  }

private:
  double chi3, fraction, gamma, lambda, R0;
};
