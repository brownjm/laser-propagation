#ifndef RADIAL_H_
#define RADIAL_H_

#include <complex>
#include <vector>
#include <fftw3.h>
#include "array2d.h"

class Radial {
public:
  Radial(int Ntime, double time_min, double time_max,
         double wave_min, double wave_max,
	 int Nradius, double Rmax, int Nkperp);
  ~Radial();

  int Ntime, Nradius, Nomega, Nkperp, index_minimum_frequency;
  double time_min, time_max, dt, wavelength_min, wavelength_max, Rmax;
  std::vector<double> time, omega, wavelength, radius, kperp;
  Array2D<std::complex<double>> temporal, spectral;
  std::vector<double> spectral_filter, temporal_filter;

  void initialize_temporal_filter(double time_filter_min, double time_filter_max);
  void initialize_spectral_filter(double wave_filter_min, double wave_filter_max);
  
  void transform_to_spectral();
  void transform_to_temporal();

  void forward_hankel();
  void forward_fft();
  void backward_fft();
  void backward_hankel();
  
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
  
  

  fftw_plan forward_plan, backward_plan;
  std::vector<std::complex<double>> Aux_fft;
  Array2D<std::complex<double>> Aux_hankel;
  Array2D<double> dht;


  void initialize_temporal_domain();
  void initialize_spectral_domain();
  void initialize_radial_domain();
  void apply_temporal_filter();
  void apply_spectral_filter();
};

#endif // RADIAL_H_
