#include "nonlinear_response.h"

class Kerr : public NonlinearResponse {
public:
  Kerr(double n2, double n0, double pressure) {
    chi3 = 4.0/3.0 * Constants::epsilon_0 * Constants::c * n2 * pressure * std::pow(n0, 2);
  }

  void calculate_response(const std::vector<double>& radius,
                          const std::vector<double>& time,
                          const Array2D<std::complex<double>>& electric_field,
                          const Array2D<double>&,
                          Array2D<std::complex<double>>& response) override {
    for (std::size_t i = 0; i < radius.size(); ++i) {
      for (std::size_t j = 0; j < time.size(); ++j) {
        double Ereal = electric_field(i, j).real();
        response(i, j) += Constants::epsilon_0 * chi3 * std::pow(Ereal, 3);
      }
    }
  }

private:
  double chi3;
};
