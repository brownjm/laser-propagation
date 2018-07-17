#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <complex>
#include <vector>
#include <fftw3.h>

class Domain {
public:
  Domain(int Ntime, double Tmax, double wave_min, double wave_max, int Nradius, double Rmax);
  ~Domain();
  std::vector<double> time, omega, radius, kperp;
  std::vector<std::complex<double>> temporal, spectral;
  int Ntime, Nradius, Nomega, Nkperp, start;
  double wave_min, wave_max;

  void transform_to_spectral();
  void transform_to_temporal();


  // convenient access functions
  inline std::complex<double> rt(int i, int j) const {
    return temporal[i*Ntime + j];
  }

  inline std::complex<double>& rt(int i, int j) {
    return temporal[i*Ntime + j];
  }

  inline std::complex<double> kw(int i, int j) const {
    return spectral[i*Nomega + j];
  }

  inline std::complex<double>& kw(int i, int j) {
    return spectral[i*Nomega + j];
  }
  
  
private:
  fftw_plan forward_plan, backward_plan;
  std::vector<std::complex<double>> Aux_fft, Aux_hankel;
  std::vector<double> dht;
};

#endif // DOMAIN_H_
