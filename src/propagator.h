#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "field.h"
#include "radial.h"
#include "array2d.h"
#include "ionization.h"
#include "simulation_data.h"
class Linear;


class Propagator {
public:
  Propagator(int Ntime, double T, double wave_min, double wave_max,
	     double filter_min, double filter_max,
	     int Nradius, double R, int Nkperp,
	     double abs_err, double rel_error, double first_step);
  ~Propagator();
  void initialize_linear(const Linear& linear, double wave0);
  void initialize_field(const Field::Field& field);
  void initialize_kerr(double n2);
  void initialize_pressure(double pressure);
  void initialize_rate(const std::string& filename, double fraction,
		       double scaling);

  void linear_step(Radial& radial, double dz);
  void linear_step(const std::complex<double>* A, Radial& radial, double dz);

  void calculate_electron_density();
  void nonlinear_step(double& z, double zi);
  void calculate_rhs(double z, const std::complex<double>* A, std::complex<double>* dA);

  SimulationData get_data();
  
  Radial field;
  Radial Pkerr;
  Radial Jnon;
  Radial Jplasma;
  Array2D<double> Rho;
  Radial RhoE;


  std::unique_ptr<Ionization::Rate> ionization_rate;
  std::vector<double> intensity, rate;
  
  int Ntime, Nradius, Nomega, Nkperp;
  double vg, n2, fraction, scaling, pressure;
  Array2D<std::complex<double>> kz, coef, A;

  // ode solver
  double z, abserr, relerr, first_step;
  gsl_odeiv2_system system;
  gsl_odeiv2_driver* driver;
};

int RHSfunction(double z, const double y[], double dy[], void* p);


#endif // PROPAGATOR_H_
