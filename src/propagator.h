#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include <sstream>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "field.h"
#include "radial.h"
#include "array2d.h"
#include "ionization.h"
#include "simulation_data.h"
#include "nonlinear_response.h"

class Linear;


class Propagator {
public:
  Propagator(int Ntime, double time_min, double time_max,
             double wave_min, double wave_max,
	     int Nradius, double R, int Nkperp,
	     double abs_err, double rel_error, double first_step);
  ~Propagator();
  std::string log_grid_info();
  void initialize_linear(const Linear& linear, double wave0);
  void initialize_field(const Field::Field& field);

  void initialize_filters(double time_filter_min, double time_filter_max,
                          double wave_filter_min, double wave_filter_max);

  void add_polarization(std::unique_ptr<NonlinearResponse> P);
  void add_current(std::unique_ptr<NonlinearResponse> J);
  void add_ionization(std::unique_ptr<Ionization::IonizationModel> ioniz);
  
  void linear_step(Radial& radial, double dz);
  void linear_step(const std::complex<double>* A, Radial& radial, double dz);

  void calculate_electron_density();
  void nonlinear_step(double& z, double zi);
  void calculate_rhs(double z, const std::complex<double>* A, std::complex<double>* dA);

  SimulationData get_data();
  
  Radial field;
  Array2D<double> Rho;

  
  int Ntime, Nradius, Nomega, Nkperp;
  double vg, n2, fraction, scaling, pressure;
  Array2D<std::complex<double>> kz, coef, A;

  std::vector<std::unique_ptr<NonlinearResponse>> polarization_responses;
  std::vector<std::unique_ptr<NonlinearResponse>> current_responses;
  std::vector<std::unique_ptr<Radial>> polarization_workspaces;
  std::vector<std::unique_ptr<Radial>> current_workspaces;
  std::unique_ptr<Ionization::IonizationModel> ionization;
  
  // ode solver
  double z, abserr, relerr, first_step;
  gsl_odeiv2_system system;
  gsl_odeiv2_driver* driver;
};

int RHSfunction(double z, const double y[], double dy[], void* p);


#endif // PROPAGATOR_H_
