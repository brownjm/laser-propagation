#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include "domain.h"
#include "constants.h"
#include "io.h"

Domain::Domain(int Nt, double Tmax, double wave_min, double wave_max, int Nr,
	       double Rmax)
  : Ntime(Nt), Nradius(Nr), Nkperp(Nr), start(0) {
  
  // construct temporal coordinates
  double dt = Tmax / (Ntime-1);
  for (int i = 0; i < Ntime; ++i) {
    double t = (i - Ntime/2) * dt;
    time.push_back(t);
  }
  temporal.resize(Ntime * Nradius);
  Aux_fft.resize(Ntime);

  forward_plan = fftw_plan_dft_1d(Ntime,
				  reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				  reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				  FFTW_FORWARD, FFTW_ESTIMATE);
  backward_plan = fftw_plan_dft_1d(Ntime,
				   reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				   reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				   FFTW_BACKWARD, FFTW_ESTIMATE);
  
  // find active frequencies
  double omega_min = 2*Constants::pi*Constants::c / wave_max;
  double omega_max = 2*Constants::pi*Constants::c / wave_min;

  for (int i = 0; i < Ntime/2; ++i) {
    double o = 2*Constants::pi * i / Tmax;
    if (o < omega_min) ++start;
    if ((o >= omega_min) && (o <= omega_max)) omega.push_back(o);
  }
  if (omega.size() < 2) {
    std::ostringstream os;
    os << "Temporal array sampling does not support requested wavelength range:\n(";
    os << std::scientific << wave_min << ", " << wave_max << ")";
    throw std::runtime_error(os.str());
  }
  this->wave_max = 2*Constants::pi*Constants::c / omega.front();
  this->wave_min = 2*Constants::pi*Constants::c / omega.back();

  Nomega = omega.size();
  Aux_hankel.resize(Nomega*Nradius);
  spectral.resize(Nomega*Nkperp);

  
  // set up radial coordinates
  std::vector<double> zeros;
  for (int i = 0; i <= Nradius; ++i) {
    zeros.push_back(gsl_sf_bessel_zero_J0(i+1));
  }

  for (int i = 0; i < Nradius; ++i) {
    radius.push_back(Rmax * zeros[i] / zeros[Nradius]);
  }
  for (int i = 0; i < Nkperp; ++i) {
    kperp.push_back(zeros[i] / Rmax);
  }

  // create transform matrix
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Nradius; ++j) {
      double J0 = gsl_sf_bessel_J0(zeros[i]*zeros[j] / zeros[Nradius]);
      double J1 = gsl_sf_bessel_J1(zeros[j]);
      double S = zeros[Nradius];
      dht.push_back(2.0 * J0 / (J1*J1*S));
    }
  }
  IO::write("hankel.dat", dht);
}


Domain::~Domain() {
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(backward_plan);
  fftw_cleanup();
}


void Domain::transform_to_spectral() {
  // Perform Fourier transform on each temporal row
  // Copy active frequencies from Aux_fft -> Aux_hankel
  for (int i = 0; i < Nradius; ++i) {
    std::copy_n(std::begin(temporal) + i*Ntime, Ntime, std::begin(Aux_fft));
    fftw_execute(forward_plan);
    std::copy_n(std::begin(Aux_fft)+start, Nomega, std::begin(Aux_hankel) + i*Nomega);
  }
  
  // Perform Hankel transform on each frequency column
  // Copy data from Aux_hankel -> spectral
  for (int o = 0; o < Nomega; ++o) {
    for (int i = 0; i < Nkperp; ++i) {
      spectral[i*Nomega + o] = 0.0;
      for (int j = 0; j < Nradius; ++j) {
	spectral[i*Nomega + o] += dht[Nradius*i + j] * Aux_hankel[j*Nomega + o];
      }
    }
  }
}


void Domain::transform_to_temporal() {
  // Hankel transform on each frequency column
  // Copy data from spectral -> Aux_hankel
  for (int o = 0; o < Nomega; ++o) {
    for (int i = 0; i < Nradius; ++i) {
      Aux_hankel[i*Nomega + o] = 0.0;
      for (int j = 0; j < Nkperp; ++j) {
	Aux_hankel[i*Nomega + o] += dht[Nradius*i + j] * spectral[j*Nomega + o];
      }
    }
  }

  // Perform inverse fourier transform on each row
  // copy data from Aux_fft -> temporal
  for (int i = 0; i < Nradius; ++i) {
    std::fill(std::begin(Aux_fft), std::end(Aux_fft), 0);
    std::copy_n(std::begin(Aux_hankel)+i*Nomega, Nomega, std::begin(Aux_fft)+start);
    fftw_execute(backward_plan);
    for (auto& a : Aux_fft) { a /= Ntime; }
    std::copy_n(std::begin(Aux_fft), Ntime, std::begin(temporal)+i*Ntime);
  }
}
