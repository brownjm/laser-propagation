#include "argon_response_pool.h"
#include "argon.h"
#include "../core/radial.h"
#include "../util/io.h"
#include "../util/constants.h"
#include <mpi.h>

ArgonResponsePool::ArgonResponsePool(int Nr, int Nl, int Nmask,
                                   const std::string& filename_potentials,
                                   double ionization_box_size,
                                   int Nradius, int Nt, double atomic_dt,
                                   double density_of_neutrals)
  :argon(Nr, Nl, Nmask, filename_potentials, ionization_box_size),
   Nradius(Nradius), Nt(Nt), atomic_dt(atomic_dt), density_of_neutrals(density_of_neutrals),
   field_atomic(Nradius, Nt), dipole(Nradius, Nt), probability_free(Nradius, Nt),
   dipole_linear(Nradius, Nt), local_field_atomic(Nt), local_dipole(Nt),
   local_probability_free(Nt), local_dipole_linear(Nt),
   temporal_filter(Nt) {
  int steps = 3000;
  double dt = 0.1;
  double loss = 1;
  argon.find_ground_state(steps, dt, loss);

  // initialize the temporal filter to ones
  std::fill(std::begin(temporal_filter), std::end(temporal_filter), 1.0);
}

void ArgonResponsePool::calculate_electron_density(const Radial& electric_field,
                                                  Array2D<double>& ionization_rate,
                                                  Array2D<double>& electron_density) {

  prepare_workers_for_electron_density();
  // both ionization and dipole moment are calculated in order to save
  // computational time
  double field_dt = electric_field.time[1] - electric_field.time[0];
  double field_atomic_dt = field_dt / au_time;

  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Nt; ++j) {
      field_atomic(i, j) = electric_field.rt(i, j).real() / au_efield;
    }
  }
  
  // send flag to keep workers alive and send them the electric field
  MPI_Scatter(field_atomic.get_data_ptr(), Nt, MPI_DOUBLE,
              local_field_atomic.data(), Nt, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // compute local dipole and free electron probability
  argon.reset_to_ground_state();
  argon.calculate_dipole_ionization(local_field_atomic, field_atomic_dt, atomic_dt,
                                    local_dipole, local_probability_free);

  // gather dipole and probability back to ArgonResponsePool
  MPI_Gather(local_dipole.data(), Nt, MPI_DOUBLE,
             dipole.get_data_ptr(), Nt, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gather(local_probability_free.data(), Nt, MPI_DOUBLE,
             probability_free.get_data_ptr(), Nt, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // force free electron probability to be constant or monotonically increasing
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 1; j < Nt; ++j) {
      if (probability_free(i, j) < probability_free(i, j-1)) {
        probability_free(i, j) = probability_free(i, j-1);
      }
    }
  }

  // set values for electron density in SI units
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Nt; ++j) {
      // argon gnd state has two electrons in m=0
      // use this prob_free instead of SAE prob.
      double prob_free = 1 - std::pow(1 - probability_free(i, j), 2);
      electron_density(i, j) = density_of_neutrals * prob_free;
    }
  }

  // calculate the rate of ionization by differentiating the free
  // electron probability using the midpoint rule
  for (int i = 0; i < Nradius; ++i) {
    ionization_rate(i, 0) = (probability_free(i, 1) - 0.0) / (2*field_dt);
    for (int j = 1; j < Nt-1; ++j) {
      ionization_rate(i, j) = (probability_free(i, j+1) - probability_free(i, j-1)) / (2*field_dt);
    }
    ionization_rate(i, Nt-1) = (0.0 - probability_free(i, Nt-2)) / (2*field_dt);
  }
}

void ArgonResponsePool::calculate_response(const std::vector<double>&,
                                          const std::vector<double>& time,
                                          const Array2D<std::complex<double>>& electric_field,
                                          const Array2D<double>&,
                                          Array2D<std::complex<double>>& response) {
  prepare_workers_for_response();
  double field_dt = time[1] - time[0];
  double field_atomic_dt = field_dt / au_time;

  // calculate the dipole moment with reduced field strength to obtain
  // linear part of dipole moment
  double scale_factor = 1e3;
  for (int j = 0; j < Nt; ++j) {
    local_field_atomic[j] = electric_field(0, j).real() / (au_efield * scale_factor);
  }

  // local computation of linear dipole moment
  argon.reset_to_ground_state();
  argon.calculate_dipole(local_field_atomic, field_atomic_dt, atomic_dt, local_dipole_linear);
  for (auto& d: local_dipole_linear) {
    d *= scale_factor;
  }
  MPI_Gather(local_dipole_linear.data(), Nt, MPI_DOUBLE,
             dipole_linear.get_data_ptr(), Nt, MPI_DOUBLE,
             0, MPI_COMM_WORLD);
             
             
  // compute nonlinear response
  for (int i = 0; i < Nradius; ++i) {
    for (int j = 0; j < Nt; ++j) {
      double nonlinear_dipole = dipole(i, j) - dipole_linear(i, j);

      // filter the dipole
      nonlinear_dipole *= temporal_filter[j];
      
      // there are two m=0 electrons in argon, multiply response by 2
      response(i, j) += 2 * density_of_neutrals * au_dipole * nonlinear_dipole;
    }
  }
}

void ArgonResponsePool::initialize_temporal_filter(std::vector<double>& time,
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



void prepare_workers_for_electron_density() {
  int flag = 1;
  MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void prepare_workers_for_response() {
  int flag = 2;
  MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void stop_workers() {
  int flag = 0;
  MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
