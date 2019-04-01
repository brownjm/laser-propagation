#include "linear.h"

namespace Linear {
  Base::Base(Medium::IndexFunction linear_index)
    :n(linear_index) {}
  
  double Base::group_velocity(double kperp, double omega) const {
    double domega = 1e-2 * omega;
    // five-point stencil to compute the derivative d/dw k(w)| w = w0
    double kp = (-kz(kperp, omega+2*domega) + 8.0*kz(kperp, omega+domega) - 8.0*kz(kperp, omega-domega) + kz(kperp, omega-2*domega)).real() / (12*domega);
    return 1/kp;
    
  }

  double Base::gvd(double kperp, double omega) const {
    double domega = 1e-2 * omega;
    // five-point stencil to compute the derivative d^2/dw^2 k(w)| w=w0
    double kpp = (-kz(kperp, omega+2*domega) + 16.0*kz(kperp, omega+domega) - 30.0*kz(kperp, omega) + 16.0*kz(kperp, omega-domega) - kz(kperp, omega-2*domega)).real() / (12*domega*domega);
    return kpp;
  }

  
  std::complex<double> FreeSpace::kz(double kperp, double omega) const {
    std::complex<double> index = n(omega);
    double k0 = omega / Constants::c;
    double k = std::real(index) * k0;
    double alpha = 2*std::imag(index) * k0;
    
    double k2 = std::pow(k, 2);
    double kperp2 = std::pow(kperp, 2);
    if (kperp2 > k2) {
      return 0.0;
    }
    else {
      double kz_real = std::sqrt(k2 - kperp2);
      std::complex<double> kzvalue(kz_real, alpha);
      return kzvalue;
    }
  }

  
  std::complex<double> DiffractionLess::kz(double, double omega) const {
    std::complex<double> index = n(omega);
    double k0 = omega / Constants::c;
    double k = std::real(index) * k0;
    double alpha = 2*std::imag(index) * k0;
    
    std::complex<double> kzvalue(k, alpha);
    return kzvalue;
  }


  std::complex<double> Capillary::kz(double kperp, double omega) const {
    double k0 = omega / Constants::c;
    std::complex<double> np = Medium::pressurize(pressure, n, omega);
    std::complex<double> k = np * k0;
    std::complex<double> beta = std::sqrt(std::pow(k, 2) - std::pow(kperp, 2));
    double eps = std::pow(nclad, 2); 
    double alpha = 1.0/(2*R) * std::pow(kperp / k0, 2) * (eps + 1) / std::sqrt(eps - 1);
    return beta + std::complex<double>(0, alpha);
  }

} // end namespace Linear
