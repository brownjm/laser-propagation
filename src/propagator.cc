#include "propagator.h"
#include "constants.h"
#include "field.h"
#include "radial.h"
#include "util.h"
#include "linear.h"
#include "io.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

Propagator::Propagator(int Nt, double time_min, double time_max,
                       double wave_min, double wave_max,
		       int Nr, double R, int Nk,
		       double abs_err, double rel_err, double step, double z)
  :field(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   electron_density(Nr, Nt), ionization_rate(Nr, Nt),
   Ntime(Nt), Nradius(Nr), vg(0),
   kz(field.Nkperp, field.Nomega),
   coef(field.Nkperp, field.Nomega),
   A(field.Nkperp, field.Nomega),
   workspace1(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   workspace2(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   current_distance(z), step(step) {


  Nomega = field.Nomega;
  Nkperp = field.Nkperp;

  std::vector<double> wavelengths;
  for (auto o : field.omega) wavelengths.push_back(2*Constants::pi*Constants::c / o);

  // set up ode solver
  system = {RHSfunction, nullptr, 2*A.vec().size(), this};
  stepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 2*A.vec().size());
  control = gsl_odeiv2_control_y_new(abs_err, rel_err);
  evolve = gsl_odeiv2_evolve_alloc(2*A.vec().size());
}

Propagator::~Propagator() {
  gsl_odeiv2_evolve_free(evolve);
  gsl_odeiv2_control_free(control);
  gsl_odeiv2_step_free(stepper);
}

std::string Propagator::log_grid_info() {
  // log message about computation box size and parameters
  double wave_min = 2*Constants::pi*Constants::c / field.omega.back();
  double wave_max = 2*Constants::pi*Constants::c / field.omega.front();
  std::stringstream ss;
  ss << "*** Computational Grid ***\n";
  ss << "Supported Wavelengths: (" << wave_min << ", " << wave_max << ")\n";
  ss << "Ntime    = " << Ntime << "\n";
  ss << "Nomega   = " << Nomega << "\n";
  ss << "Nradius  = " << Nradius << "\n";
  ss << "Nkperp   = " << Nkperp << "\n";

  // ode solver
  ss << "ODE size = " << A.vec().size() << " complex<double>\n\n";
  return ss.str();
}

void Propagator::initialize_linear(const Linear::Base& linear, double omega0) {
  for (int j = 0; j < Nomega; ++j) {
    index.push_back(linear.n(field.omega[j]));
  }
  
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    double kperp = field.kperp[i];
    for (int j = 0; j < Nomega; ++j) {
      double omega = field.omega[j];
      kz(i, j) = linear.kz(kperp, omega);
      if (std::isnan(kz(i, j).real())) {
        std::cout << omega << " " << kperp << "\n";
      }
    }
  }
  vg = linear.group_velocity(field.kperp[0], omega0);

  // nonlinear coupling coefficient
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      double kzvalue = kz(i, j).real();
      if (kzvalue > 0.0) {
        coef(i, j) = field.omega[j] / (2*Constants::epsilon_0*std::pow(Constants::c, 2)*kz(i, j));
      }
      else {
        coef(i, j) = 0.0;
      }
    }
  }
}


void Propagator::initialize_field(const Field::Field& Efield) {
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Ntime; ++j) {
      field.rt(i, j) = Efield(current_distance, field.radius[i], field.time[j]);
    }
  }
  
  // fill spectral array with initialized field
  field.transform_to_spectral();

  // copy spectral field to auxillary A which is passed to ODE solver
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      A(i, j) = field.ko(i, j);
    }
  }

  field.transform_to_temporal();
}

void Propagator::restart_from(const std::string& spectral_filename) {
  std::vector<std::complex<double>> spectral;
  IO::read_binary(spectral_filename, spectral);
  if (spectral.size() == field.spectral.values.size()) {
    field.spectral.values = spectral;
  }
  else {
    throw std::runtime_error("Spectral field file dimensions do not match the current simulation parameters. Expected a length of " + std::to_string(Nkperp*Nomega) + ", received " + std::to_string(spectral.size()));
  }
  
  // copy spectral field to auxillary A which is passed to ODE solver
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      A(i, j) = field.ko(i, j);
    }
  }

  field.transform_to_temporal();
}

void Propagator::initialize_filters(double time_filter_min, double time_filter_max,
                                    double wave_filter_min, double wave_filter_max) {
  field.initialize_temporal_filter(time_filter_min, time_filter_max);
  field.initialize_spectral_filter(wave_filter_min, wave_filter_max);
  workspace1.initialize_temporal_filter(time_filter_min, time_filter_max);
  workspace1.initialize_spectral_filter(wave_filter_min, wave_filter_max);
  workspace2.initialize_temporal_filter(time_filter_min, time_filter_max);
  workspace2.initialize_spectral_filter(wave_filter_min, wave_filter_max);
}


void Propagator::add_polarization(std::shared_ptr<NonlinearResponse> polarization) {
  polarization_responses.push_back(polarization);
}

void Propagator::add_current(std::shared_ptr<NonlinearResponse> current) {
  current_responses.push_back(current);
}

void Propagator::add_ionization(std::shared_ptr<Ionization> ioniz) {
  ionization = ioniz;
}

void Propagator::linear_step(Radial& radial, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      auto arg = kz(i, j) - radial.omega[j] / vg;
      radial.ko(i, j) *= std::exp(imagi * arg * dz);
    }
  }
}

void Propagator::linear_step(std::complex<double>* A, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      auto arg = kz(i, j) - field.omega[j] / vg;
      A[i*Nomega + j] *= std::exp(imagi * arg * dz);
    }
  }
}

void Propagator::linear_step(const std::complex<double>* A, Radial& radial, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      auto arg = kz(i, j) - radial.omega[j] / vg;
      radial.ko(i, j) = A[i*Nomega + j] * std::exp(imagi * arg * dz);
    }
  }
}

void Propagator::nonlinear_step(double& z, double z_end) {
  current_distance = z;
  double last_step = 0;

  while (z < z_end) {
    int status = gsl_odeiv2_evolve_apply(evolve, control, stepper,
                                         &system,
                                         &z, z_end,
                                         &step,
                                         reinterpret_cast<double*>(A.get_data_ptr()));

    if (status != GSL_SUCCESS) {
      throw std::runtime_error("gsl_ode error: " + std::to_string(status));
    }

    last_step = z - current_distance;
    current_distance = z;

    // linearly propagate spectral field
    linear_step(A.get_data_ptr(), last_step);
  }

  // copy solver's A to the field used for calculating observables
  field.spectral.values = A.values;
  field.transform_to_temporal();
  calculate_electron_density();
}

void Propagator::calculate_electron_density() {
  if (ionization) {
    ionization->calculate_electron_density(field, ionization_rate, electron_density);
  }
}

void Propagator::calculate_rhs(double z, const std::complex<double>* A, std::complex<double>* dA) {
  // proposed step size
  double dz = z - current_distance;
  
  // apply linear propagation over distance dz
  linear_step(A, field, dz);

  // transform field(k, omega) -> field(r, t) and calculate electron density rho(r, t)
  field.transform_to_temporal();
  calculate_electron_density();

  // calculate contributions from nonlinear polarizations Pnl(r, t)
  std::fill(std::begin(workspace1.temporal.values), std::end(workspace1.temporal.values), 0);
  for (auto& source : polarization_responses) {
    source->calculate_response(field.radius, field.time, field.temporal, electron_density,
                               workspace1.temporal);
  }
  // transform Pnl(r, t) -> Pnl(r, omega)
  workspace1.backward_fft();

  // calculate contributions from nonlinear currents Jnl(r, t)
  std::fill(std::begin(workspace2.temporal.values), std::end(workspace2.temporal.values), 0);
  for (auto& source : current_responses) {
    source->calculate_response(field.radius, field.time, field.temporal, electron_density,
                               workspace2.temporal);
  }
  // transform Jnl(r, t) -> Jnl(r, omega)
  workspace2.backward_fft();

  // combine responses: i omega Pnl(r, omega) - Jnl(r, omega)
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      double omega = field.omega[j];
      workspace1.ro(i, j) = imagi*omega*workspace1.ro(i, j) - workspace2.ro(i, j);
    }
  }

  // transform total response to (k, omega)
  workspace1.backward_hankel();

  // multiply total response (k, omega) by nonlinear coupling
  // coefficient and place into dA
  std::transform(std::begin(workspace1.spectral.values),
                 std::end(workspace1.spectral.values),
                 std::begin(coef.values), dA, std::multiplies<std::complex<double>>());

  // apply linear propagation to undo previous shift dz
  linear_step(dA, -dz);
}


SimulationData Propagator::get_data() {
  SimulationData data = {field, electron_density, *this};
  return data;
}

int RHSfunction(double z, const double y[], double dy[], void* p) {
  Propagator& prop = *static_cast<Propagator*>(p);
  const std::complex<double>* A = reinterpret_cast<const std::complex<double>*>(y);
  std::complex<double>* dA = reinterpret_cast<std::complex<double>*>(dy);
  prop.calculate_rhs(z, A, dA);
  return GSL_SUCCESS;
}

