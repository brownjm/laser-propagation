#ifndef RADIAL_H_
#define RADIAL_H_

#include <complex>
#include <vector>
#include <fftw3.h>
#include "array2d.h"

class Radial {
public:
  Radial(int Ntime, double Tmax, double wave_min, double wave_max,
	 double filter_min, double filter_max,
	 int Nradius, double Rmax, int Nkperp);
  ~Radial();

  int Ntime, Nradius, Nomega, Nkperp, index_minimum_frequency;
  double Tmax, dt, wavelength_min, wavelength_max, Rmax;
  std::vector<double> time, omega, wavelength, radius, kperp;
  Array2D<std::complex<double>> temporal, spectral;
  std::vector<double> spectral_filter, temporal_filter;
  
  void transform_to_spectral();
  void transform_to_temporal();


  // convenient access functions
  inline std::complex<double> rt(int i, int j) const {
    return temporal(i, j);
  }

  inline std::complex<double>& rt(int i, int j) {
    return temporal(i, j);
  }

  inline std::complex<double> kw(int i, int j) const {
    return spectral(i, j);
  }

  inline std::complex<double>& kw(int i, int j) {
    return spectral(i, j);
  }
  
  
private:
  fftw_plan forward_plan, backward_plan;
  std::vector<std::complex<double>> Aux_fft;
  Array2D<std::complex<double>> Aux_hankel;
  Array2D<double> dht;


  void initialize_temporal_domain();
  void initialize_temporal_filter();
  void initialize_spectral_domain();
  void initialize_spectral_filter(double filter_min, double filter_max);
  void initialize_radial_domain();
  void apply_temporal_filter();
  void apply_spectral_filter();
};

#endif // RADIAL_H_
