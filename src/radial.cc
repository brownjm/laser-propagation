#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include "radial.h"
#include "constants.h"
#include "io.h"
#include <iostream>



Radial::Radial(int Nt, double T, double wave_min, double wave_max,
	       double filter_min, double filter_max,
	       int Nr, double R, int Nk)
  : Ntime(Nt), Nradius(Nr), Nomega(0), Nkperp(Nk), index_minimum_frequency(0),
    Tmax(T), dt(Tmax / (Ntime-1)), wavelength_min(wave_min), wavelength_max(wave_max), Rmax(R) {

  if (Nkperp > Nradius) {
    std::ostringstream os;
    os << "Nkperp (" << Nkperp << ") cannot be greater than Nradius (" << Nradius << ")";
    throw std::runtime_error(os.str());
  }

  initialize_temporal_domain();
  initialize_temporal_filter();
  initialize_spectral_domain();
  initialize_spectral_filter(filter_min, filter_max);
  initialize_radial_domain();
}

void Radial::initialize_temporal_domain() {
  // construct temporal coordinates
  for (int i = 0; i < Ntime; ++i) {
    double t = (i - Ntime/2) * dt;
    time.push_back(t);
  }
  temporal.resize(Nradius, Ntime);
  Aux_fft.resize(Ntime);

  forward_plan = fftw_plan_dft_1d(Ntime,
				  reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				  reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				  FFTW_FORWARD, FFTW_MEASURE);
  //FFTW_FORWARD, FFTW_ESTIMATE);
  backward_plan = fftw_plan_dft_1d(Ntime,
				   reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				   reinterpret_cast<fftw_complex*>(Aux_fft.data()),
				   FFTW_BACKWARD, FFTW_MEASURE);
  //FFTW_BACKWARD, FFTW_ESTIMATE);

  // fftw plan may have written data to the auxillary array if
  // FFTW_{MEASURE,PATIENT,EXHAUSTIVE} is used, therefore we set all values to zero
  std::fill(std::begin(Aux_fft), std::end(Aux_fft), 0);
}

void Radial::initialize_spectral_domain() {
  // find active frequencies
  double omega_min = 2*Constants::pi*Constants::c / wavelength_max;
  double omega_max = 2*Constants::pi*Constants::c / wavelength_min;

  index_minimum_frequency = 0;
  for (int i = 0; i < Ntime/2; ++i) {
    double o = 2*Constants::pi * i / Tmax;
    if (o < omega_min) ++index_minimum_frequency;
    if ((o >= omega_min) && (o <= omega_max)) omega.push_back(o);
  }
  if (omega.size() < 2) {
    std::ostringstream os;
    os << "Temporal array sampling dt = ";
    os << dt;
    os << " does not support requested wavelength range:\n(";
    os << std::scientific << wavelength_min << ", " << wavelength_max << ")";
    throw std::runtime_error(os.str());
  }

  for (auto o : omega) wavelength.push_back(2*Constants::pi*Constants::c / o);
  this->wavelength_max = wavelength.front();
  this->wavelength_min = wavelength.back();

  Nomega = omega.size();
  Aux_hankel.resize(Nradius, Nomega);
  spectral.resize(Nkperp, Nomega);
}

void Radial::initialize_spectral_filter(double filter_min, double filter_max) {
  // set up filters
  spectral_filter.resize(Nomega);
  std::fill(std::begin(spectral_filter), std::end(spectral_filter), 1.0);

  double filter_omega_min = 2*Constants::pi*Constants::c / filter_max;
  double filter_omega_max = 2*Constants::pi*Constants::c / filter_min;
  int low = 0;
  int high = 0;
  for (int j = 0; j < Nomega; ++j) {
    if (omega[j] < filter_omega_min) ++low;
    if (omega[j] < filter_omega_max) ++high;
  }
  // std::cout << "low  = " << low << "\n";
  // std::cout << "high = " << high << "\n";

  for (int j = 0; j < low; ++j) {
    spectral_filter[j] = std::pow(std::sin(Constants::pi/2 * j / (low-1)), 2);
  }

  for (int j = high; j < Nomega; ++j) {
    spectral_filter[j] = std::pow(std::cos(Constants::pi/2 * (j - high) / (Nomega-high-1)), 2);
  }
  
  IO::write("spectral_filter.dat", spectral_filter);
}

void Radial::initialize_temporal_filter() {
  temporal_filter.resize(Ntime);
  std::fill(std::begin(temporal_filter), std::end(temporal_filter), 1.0);
  int width = Ntime / 5;
  for (int j = 0; j < width; ++j) {
    double w = std::pow(std::sin(Constants::pi/2 * j / (width-1)), 2);
    temporal_filter[j] = w;
    temporal_filter[Ntime-1 - j] = w;
  }

  IO::write("temporal_filter.dat", temporal_filter);
}

void Radial::initialize_radial_domain() {
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
  dht.resize(Nradius, Nradius);
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Nradius; ++j) {
      double J0 = gsl_sf_bessel_J0(zeros[i]*zeros[j] / zeros[Nradius]);
      double J1 = gsl_sf_bessel_J1(zeros[j]);
      double S = zeros[Nradius];
      dht(i, j) = 2.0 * J0 / (J1*J1*S);
    }
  }
  //IO::write("hankel.dat", dht.vec());
}


Radial::~Radial() {
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(backward_plan);
  fftw_cleanup();
}


void Radial::apply_temporal_filter() {
  // apply temporal guard
  for (int j = 0; j < Ntime; ++j) {
    Aux_fft[j] *= temporal_filter[j];
  }
}

void Radial::apply_spectral_filter() {
  // apply spectral filter
  for (int j = 0; j < Nomega; ++j) {
    Aux_fft[j+index_minimum_frequency] *= spectral_filter[j];
  }
}

void Radial::transform_to_spectral() {
  // Perform Fourier transform on each temporal row
  // Copy active frequencies from Aux_fft -> Aux_hankel
  for (int i = 0; i < Nradius; ++i) {
    std::fill(std::begin(Aux_fft), std::end(Aux_fft), 0);
    std::copy_n(std::begin(temporal.values) + i*Ntime, Ntime, std::begin(Aux_fft));
    apply_temporal_filter();    
    fftw_execute(forward_plan);
    apply_spectral_filter();
    std::copy_n(std::begin(Aux_fft)+index_minimum_frequency,
		Nomega, std::begin(Aux_hankel.values) + i*Nomega);
  }
  
  // Perform Hankel transform on each frequency column
  // Copy data from Aux_hankel -> spectral
  for (int o = 0; o < Nomega; ++o) {
    for (int i = 0; i < Nkperp; ++i) {
      spectral(i, o) = 0.0;
      for (int j = 0; j < Nradius; ++j) {
	spectral(i, o) += dht(i, j) * Aux_hankel(j, o);
      }
    }
  }
}


void Radial::transform_to_temporal() {
  // Hankel transform on each frequency column
  // Copy data from spectral -> Aux_hankel
  for (int o = 0; o < Nomega; ++o) {
    for (int i = 0; i < Nradius; ++i) {
      Aux_hankel(i, o) = 0.0;
      for (int j = 0; j < Nkperp; ++j) {
	Aux_hankel(i, o) += dht(i, j) * spectral(j, o);
      }
    }
  }

  // Perform inverse fourier transform on each row
  // copy data from Aux_fft -> temporal
  for (int i = 0; i < Nradius; ++i) {
    std::fill(std::begin(Aux_fft), std::end(Aux_fft), 0);
    std::copy_n(std::begin(Aux_hankel.values) + i*Nomega, Nomega,
		std::begin(Aux_fft) + index_minimum_frequency);
    apply_spectral_filter();
    fftw_execute(backward_plan);
    apply_temporal_filter();
    for (auto& a : Aux_fft) { a /= Ntime; }
    std::copy_n(std::begin(Aux_fft), Ntime, std::begin(temporal.values)+i*Ntime);
  }
}



