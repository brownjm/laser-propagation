#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include "radial.h"
#include "constants.h"
#include "io.h"
#include <iostream>

#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Md;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mcd;


Radial::Radial(int Nt, double time_min, double time_max,
               double wave_min, double wave_max,
	       int Nr, double R, int Nk)
  : Ntime(Nt), Nradius(Nr), Nomega(0), Nkperp(Nk), index_minimum_frequency(0),
    time_min(time_min), time_max(time_max),
    wavelength_min(wave_min), wavelength_max(wave_max), Rmax(R) {
  dt = (time_max - time_min) / (Ntime - 1);
  if (Nkperp > Nradius) {
    std::ostringstream os;
    os << "Nkperp (" << Nkperp << ") cannot be greater than Nradius (" << Nradius << ")";
    throw std::runtime_error(os.str());
  }

  initialize_temporal_domain();
  initialize_spectral_domain();
  initialize_radial_domain();
}

void Radial::initialize_temporal_domain() {
  // construct temporal coordinates
  for (int i = 0; i < Ntime; ++i) {
    double t = i*dt + time_min;
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

  // fill temporal filter values to unity
  temporal_filter.resize(Ntime);
  std::fill(std::begin(temporal_filter), std::end(temporal_filter), 1.0);
}

void Radial::initialize_spectral_domain() {
  // find active frequencies
  double omega_min = 2*Constants::pi*Constants::c / wavelength_max;
  double omega_max = 2*Constants::pi*Constants::c / wavelength_min;

  index_minimum_frequency = 0;
  for (int i = 0; i < Ntime/2; ++i) {
    double o = 2*Constants::pi * i / (time_max - time_min);
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

  // set spectral filter values to unity
  spectral_filter.resize(Nomega);
  std::fill(std::begin(spectral_filter), std::end(spectral_filter), 1.0);
}

void Radial::initialize_spectral_filter(double wave_filter_min, double wave_filter_max) {
  double omega_filter_min = 2*Constants::pi*Constants::c / wave_filter_max;
  double omega_filter_max = 2*Constants::pi*Constants::c / wave_filter_min;
  int low = 0;
  int high = 0;
  for (auto o : omega) {
    if (o < omega_filter_min) ++low;
    if (o < omega_filter_max) ++high;
  }

  for (int j = 0; j < low; ++j) {
    double sin = std::sin(Constants::pi/2 * j / (low-1));
    spectral_filter[j] = std::pow(sin, 2);
  }

  for (int j = high; j < Nomega; ++j) {
    double cos = std::cos(Constants::pi/2 * (j - high) / (Nomega-high-1));
    spectral_filter[j] = std::pow(cos, 2);
  }
  
  IO::write("spectral_filter.dat", spectral_filter);
}

void Radial::initialize_temporal_filter(double time_min, double time_max) {
  int low = 0;
  int high = 0;
  for (auto t : time) {
    if (t < time_min) ++low;
    if (t < time_max) ++high;
  }

  for (int j = 0; j < low; ++j) {
    double sin = std::sin(Constants::pi/2 * j / (low-1));
    temporal_filter[j] = std::pow(sin, 2);
  }

  for (int j = high; j < Ntime; ++j) {
    double cos = std::cos(Constants::pi/2 * (j - high) / (Ntime-high-1));
    temporal_filter[j] = std::pow(cos, 2);
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
  IO::write("hankel.dat", dht.vec(), Nradius, Nradius);
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
  // for (int o = 0; o < Nomega; ++o) {
  //   for (int i = 0; i < Nkperp; ++i) {
  //     spectral(i, o) = 0.0;
  //     for (int j = 0; j < Nradius; ++j) {
  // 	spectral(i, o) += dht(i, j) * Aux_hankel(j, o);
  //     }
  //   }
  // }


  // Use eigen for Hankel transform
  Eigen::Map<Md> dht2(dht.get_data_ptr(), Nradius, Nradius);
  Eigen::Map<Mcd> Aux_hankel2(Aux_hankel.get_data_ptr(), Nradius, Nomega);
  Eigen::Map<Mcd> spectral2(spectral.get_data_ptr(), Nkperp, Nomega);
  spectral2 = dht2.block(0, 0, Nkperp, Nradius) * Aux_hankel2;
}


void Radial::transform_to_temporal() {
  // Hankel transform on each frequency column
  // Copy data from spectral -> Aux_hankel
  // for (int o = 0; o < Nomega; ++o) {
  //   for (int i = 0; i < Nradius; ++i) {
  //     Aux_hankel(i, o) = 0.0;
  //     for (int j = 0; j < Nkperp; ++j) {
  // 	Aux_hankel(i, o) += dht(i, j) * spectral(j, o);
  //     }
  //   }
  // }

  // Use eigen for Hankel transform
  Eigen::Map<Md> dht2(dht.get_data_ptr(), Nradius, Nradius);
  Eigen::Map<Mcd> Aux_hankel2(Aux_hankel.get_data_ptr(), Nradius, Nomega);
  Eigen::Map<Mcd> spectral2(spectral.get_data_ptr(), Nkperp, Nomega);
  Aux_hankel2 = dht2.block(0, 0, Nradius, Nkperp) * spectral2;

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
