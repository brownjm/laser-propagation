#ifndef GENERATE_RATE_H_
#define GENERATE_RATE_H_

class GenerateRate {
public:
  GenerateRate(double Ip, double wavelength);

  double adk(double intensity) const;
  double mpi(double intensity) const;
  double tunnel(double intensity) const;
  double yudin(double intensity) const;
  double ilkov(double intensity) const;

private:
  const double au_energy = 4.3597441775e-18;
  const double au_time = 2.41888432650516e-17;
  const double au_field = 5.1422065211e11;
  int Z, n, l, m;
  double omega, E0, F0, nstar, lstar;
  
  int factorial(int n) const;
  double f(double l, double m) const;
  double C(double n, double l) const;
  double g(double gamma) const;
  double beta(double gamma) const;
  double alpha(double gamma) const;
  double dawson(double x) const;
  double A(double E0, double omega, double gamma) const;
};


#endif // GENERATE_RATE_H_
