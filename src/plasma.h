#include "nonlinear_response.h"

class Plasma : public NonlinearResponse {
public:
  Plasma(double collision_time, double pressure)
    : collision_time(collision_time), pressure(pressure) {}

  void calculate_temporal_response(const Radial& electric_field, const Array2D<double>& electron_density, Radial& response) override {
    for (int i = 0; i < response.Nradius; ++i) {
      for (int j = 0; j < response.Ntime; ++j) {
        const double E = electric_field.rt(i, j).real();
        response.rt(i, j) = E * electron_density(i, j);
      }
    }
  }

  void finalize_spectral_response(Radial& response) override {
    const double tau = collision_time / pressure;
    const double A = std::pow(Constants::e, 2) * tau / Constants::m_e;
    const std::complex<double> imagi(0, 1);
    for (int i = 0; i < response.Nkperp; ++i) {
      for (int j = 0; j < response.Nomega; ++j) {
        const double ot = response.omega[j] * tau;
        const std::complex<double> B = (1.0 + imagi*ot) / (1.0 + std::pow(ot, 2));
        response.kw(i, j) *= A * B;
      }
    }
  }
};
