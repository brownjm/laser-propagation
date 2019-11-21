#include "argon_response.h"
#include "argon.h"
#include "../core/radial.h"
#include "../util/io.h"
#include "../util/constants.h"

ArgonResponse::ArgonResponse(int Nr, int Nl, int Nmask,
                             const std::string& filename_potentials,
                             double ionization_box_size,
                             int Nradius, int Nt, double atomic_dt,
                             double density_of_neutrals)
  :argon(Nr, Nl, Nmask, filename_potentials, ionization_box_size),
   Nradius(Nradius), Nt(Nt), atomic_dt(atomic_dt), density_of_neutrals(density_of_neutrals),
   field_atomic(Nt), probability_free(Nt), dipole_atomic(Nt),
   dipole(Nradius, Nt), temporal_filter(Nt) {
  int steps = 3000;
  double dt = 0.1;
  double loss = 1;
  argon.find_ground_state(steps, dt, loss);

  // initialize the temporal filter to ones
  std::fill(std::begin(temporal_filter), std::end(temporal_filter), 1.0);
}

void ArgonResponse::calculate_electron_density(const Radial& electric_field,
                                               Array2D<double>& ionization_rate,
                                               Array2D<double>& electron_density) {
  // both ionization and dipole moment are calculated in order to save
  // computational time
  double field_dt = electric_field.time[1] - electric_field.time[0];
  double field_atomic_dt = field_dt / au_time;
  
  // for each radial grid point
  for (int i = 0; i < Nradius; ++i) {
    // convert field to atomic units
    for (int j = 0; j < Nt; ++j) {
      field_atomic[j] = electric_field.rt(i, j).real() / au_efield;
    }
    argon.reset_to_ground_state();
    argon.calculate_dipole_ionization(field_atomic, field_atomic_dt, atomic_dt,
                                      dipole_atomic, probability_free);
    // force free electron probability to be constant or monotonically increasing
    for (int j = 1; j < Nt; ++j) {
      if (probability_free[j] < probability_free[j-1]) {
        probability_free[j] = probability_free[j-1];
      }
    }
    
    for (int j = 0; j < Nt; ++j) {
      // argon gnd state has two electrons in m=0
      // use this prob_free instead of SAE prob.
      double prob_free = 1 - std::pow(1 - probability_free[j], 2);
      electron_density(i, j) = density_of_neutrals * prob_free;
      dipole(i, j) = dipole_atomic[j];
    }

    // calculate rate of ionization by differentiating ionization using midpoint rule
    ionization_rate(i, 0) = 0.0; // assume rate starts at zero
    for (int j = 1; j < Nt-1; ++j) {
      ionization_rate(i, j) = (probability_free[j+1] - probability_free[j-1]) / (2*field_dt);
    }
    ionization_rate(i, Nt-1) = ionization_rate(i, Nt-2); // last value equal to previous value
  }
}

void ArgonResponse::calculate_response(const std::vector<double>& radius,
                                       const std::vector<double>& time,
                                       const Array2D<std::complex<double>>& electric_field,
                                       const Array2D<double>&,
                                       Array2D<std::complex<double>>& response) {
  double field_dt = time[1] - time[0];
  double field_atomic_dt = field_dt / au_time;

  // calculate the dipole moment with reduced field strength to obtain
  // linear part of dipole moment
  double scale_factor = 1e3;
  for (std::size_t i = 0; i < radius.size(); ++i) {
    for (std::size_t j = 0; j < time.size(); ++j) {
      field_atomic[j] = electric_field(i, j).real() / (au_efield * scale_factor);
    }
    argon.reset_to_ground_state();
    argon.calculate_dipole(field_atomic, field_atomic_dt, atomic_dt, dipole_atomic);
    for (std::size_t j = 0; j < time.size(); ++j) {
      double nonlinear_dipole = dipole(i, j) - scale_factor * dipole_atomic[j];

      // filter the dipole
      nonlinear_dipole *= temporal_filter[j];

      // there are two m=0 electrons in argon, multiply response by 2
      response(i, j) += 2 * density_of_neutrals * au_dipole * nonlinear_dipole;
    }
  }
}

void ArgonResponse::save_wavefunction(const std::string& filename) {
  argon.save_wavefunction(filename);
}

void ArgonResponse::initialize_temporal_filter(std::vector<double>& time,
                                               double filter_start_time) {
  int index = 0;
  for (auto t : time) {
    if (t < filter_start_time) ++index;
  }

  for (int j = index; j < Nt; ++j) {
    double cos = std::cos(Constants::pi/2 * (j - index) / (Nt-index-1));
    temporal_filter[j] = std::pow(cos, 2);
  }
}

