#include "argon_response.h"
#include "argon.h"
#include "radial.h"
#include "io.h"

ArgonResponse::ArgonResponse(int Nr, int Nl, int Nmask,
                             const std::string& filename_potentials,
                             double ionization_box_size,
                             int Nradius, int Nt, double atomic_dt,
                             double density_of_neutrals)
  :argon(Nr, Nl, Nmask, filename_potentials, ionization_box_size),
   Nradius(Nradius), Nt(Nt), atomic_dt(atomic_dt), density_of_neutrals(density_of_neutrals),
   field_atomic(Nt), probability_free(Nt), dipole_atomic(Nt),
   dipole(Nradius, Nt) {
  int steps = 3000;
  double dt = 0.1;
  double loss = 1;
  argon.find_ground_state(steps, dt, loss);
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
      electron_density(i, j) = density_of_neutrals * probability_free[j];
      dipole(i, j) = dipole_atomic[j];
    }
    // calculate rate of ionization by differentiating ionization using midpoint rule
    ionization_rate(i, 0) = (probability_free[1] - 0.0) / (2*field_dt);
    for (int j = 1; j < Nt-1; ++j) {
      ionization_rate(i, j) = (probability_free[j+1] - probability_free[j-1]) / (2*field_dt);
    }
    ionization_rate(i, Nt-1) = (0.0 - probability_free[Nt-2]) / (2*field_dt);
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
      response(i, j) += density_of_neutrals * au_dipole * nonlinear_dipole;
    }
  }
}
