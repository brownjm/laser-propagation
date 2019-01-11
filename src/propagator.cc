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

Propagator::Propagator(int Nt, double time_min, double time_max,
                       double wave_min, double wave_max,
		       int Nr, double R, int Nk,
		       double abs_err, double rel_err, double first_step)
  :field(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   electron_density(Nr, Nt),
   Ntime(Nt), Nradius(Nr), vg(0),
   kz(field.Nkperp, field.Nomega),
   coef(field.Nkperp, field.Nomega),
   A(field.Nkperp, field.Nomega),
   nonlinear_workspace(Nt, time_min, time_max, wave_min, wave_max, Nr, R, Nk),
   abserr(abs_err), relerr(rel_err), first_step(first_step) {


  Nomega = field.Nomega;
  Nkperp = field.Nkperp;

  std::vector<double> wavelengths;
  for (auto o : field.omega) wavelengths.push_back(2*Constants::pi*Constants::c / o);

  system = {RHSfunction, nullptr, 2*A.vec().size(), this};

  driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, first_step,
  					 abserr, relerr);
}

Propagator::~Propagator() {
  gsl_odeiv2_driver_free(driver);
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

void Propagator::initialize_linear(const Linear& linear, double omega0) {
  std::vector<double> index;
  for (int j = 0; j < Nomega; ++j) {
    index.push_back(linear.n(field.omega[j]));
  }
  IO::write("index.dat", index);
  
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
        coef(i, j) = imagi / (2*Constants::epsilon_0*kz(i, j).real()) * std::pow(field.omega[j]/Constants::c, 2);
      }
      else {
        coef(i, j) = 0.0;
      }
    }
  }

  IO::write("kz.dat", kz.vec(), Nkperp, Nomega);
  IO::write("coef.dat", coef.vec(), Nkperp, Nomega);
}


void Propagator::initialize_field(const Field::Field& Efield) {
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Ntime; ++j) {
      field.rt(i, j) = Efield(field.radius[i], field.time[j]);
    }
  }
  
  // fill spectral array with initialized field
  field.transform_to_spectral();

  // copy spectral field to auxillary A which is passed to ODE solver
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      A(i, j) = field.kw(i, j);
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
      A(i, j) = field.kw(i, j);
    }
  }

  field.transform_to_temporal();
}

void Propagator::initialize_filters(double time_filter_min, double time_filter_max,
                                    double wave_filter_min, double wave_filter_max) {
  field.initialize_temporal_filter(time_filter_min, time_filter_max);
  field.initialize_spectral_filter(wave_filter_min, wave_filter_max);
  nonlinear_workspace.initialize_temporal_filter(time_filter_min, time_filter_max);
  nonlinear_workspace.initialize_spectral_filter(wave_filter_min, wave_filter_max);
}


void Propagator::add_polarization(std::unique_ptr<NonlinearResponse> polarization) {
  polarization_responses.push_back(std::move(polarization));
}

void Propagator::add_current(std::unique_ptr<NonlinearResponse> current) {
  current_responses.push_back(std::move(current));
}

void Propagator::add_ionization(std::shared_ptr<Ionization::IonizationModel> ioniz) {
  ionization = ioniz;
}

void Propagator::linear_step(Radial& radial, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      auto arg = kz(i, j) - radial.omega[j] / vg;
      radial.kw(i, j) *= std::exp(-imagi * arg * dz);
    }
  }
}

void Propagator::linear_step(const std::complex<double>* A, Radial& radial, double dz) {
  std::complex<double> imagi(0, 1);
  for (int i = 0; i < Nkperp; ++i) {
    for (int j = 0; j < Nomega; ++j) {
      auto arg = kz(i, j) - radial.omega[j] / vg;
      radial.kw(i, j) = A[i*Nomega + j] * std::exp(-imagi * arg * dz);
    }
  }
}

void Propagator::nonlinear_step(double& z, double zi) {
  int status = gsl_odeiv2_driver_apply(driver, &z, zi,
        			       reinterpret_cast<double*>(A.get_data_ptr()));
  if (status != GSL_SUCCESS) {
    throw std::runtime_error("gsl_ode error: " + std::to_string(status));
  }
  linear_step(A.get_data_ptr(), field, z);
  field.transform_to_temporal();
  calculate_electron_density();
}

void Propagator::calculate_electron_density() {
  if (ionization) {
    ionization->calculate_electron_density(field, electron_density);
  }
}

void Propagator::calculate_rhs(double z, const std::complex<double>* A, std::complex<double>* dA) {
  std::fill(dA, dA + Nkperp*Nomega, 0);
  
  // 1: shift to current z
  linear_step(A, field, z);

  // 2: transform A to E
  field.transform_to_temporal();
  calculate_electron_density();

  // 3: calculate nonlinearities
  for (auto& source : polarization_responses) {
    source->calculate_temporal_response(field, electron_density, nonlinear_workspace);
    nonlinear_workspace.transform_to_spectral();
    linear_step(nonlinear_workspace, -z);
    source->finalize_spectral_response(nonlinear_workspace);
    for (int i = 0; i < Nkperp; ++i) {
      for (int j = 0; j < Nomega; ++j) {
        dA[i*Nomega + j] += coef(i, j) * nonlinear_workspace.kw(i, j);
      }
    }
  }


  std::complex<double> imagi(0, 1);
  for (auto& source : current_responses) {
    source->calculate_temporal_response(field, electron_density, nonlinear_workspace);
    nonlinear_workspace.transform_to_spectral();
    linear_step(nonlinear_workspace, -z);
    source->finalize_spectral_response(nonlinear_workspace);
    for (int i = 0; i < Nkperp; ++i) {
      for (int j = 0; j < Nomega; ++j) {
        const double omega = field.omega[j];
        dA[i*Nomega + j] += coef(i, j) * imagi / omega * nonlinear_workspace.kw(i, j);
      }
    }
  }
}


SimulationData Propagator::get_data() {
  SimulationData data = {field, electron_density};
  return data;
}

int RHSfunction(double z, const double y[], double dy[], void* p) {
  Propagator& prop = *static_cast<Propagator*>(p);
  const std::complex<double>* A = reinterpret_cast<const std::complex<double>*>(y);
  std::complex<double>* dA = reinterpret_cast<std::complex<double>*>(dy);
  prop.calculate_rhs(z, A, dA);
  return GSL_SUCCESS;
}

