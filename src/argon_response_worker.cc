#include "argon_response_worker.h"
#include <mpi.h>

ArgonResponseWorker::ArgonResponseWorker(int Nr, int Nl, int Nmask,
                                         const std::string& filename_potentials,
                                         double ionization_box_size,
                                         int Nt, double atomic_dt,
                                         double field_atomic_dt)
  :argon(Nr, Nl, Nmask, filename_potentials, ionization_box_size),
   Nt(Nt), atomic_dt(atomic_dt), field_atomic_dt(field_atomic_dt),
   field_atomic(Nt), dipole(Nt), probability_free(Nt), linear_dipole(Nt) {
  int steps = 3000;
  double dt = 0.1;
  double loss = 1;
  argon.find_ground_state(steps, dt, loss);
}

void ArgonResponseWorker::run() {
  int task;
  while (true) {
    // get task from ArgonResponsePool
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (task == 1) { // calculate electron density
      MPI_Scatter(nullptr, 0, MPI_DOUBLE,
                  field_atomic.data(), Nt, MPI_DOUBLE,
                  0, MPI_COMM_WORLD);

      // calculate dipole and ionization probability
      argon.reset_to_ground_state();
      argon.calculate_dipole_ionization(field_atomic, field_atomic_dt, atomic_dt,
                                        dipole, probability_free);

      // send data back to ArgonResponsePool
      MPI_Gather(dipole.data(), Nt, MPI_DOUBLE,
                 nullptr, 0, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

      MPI_Gather(probability_free.data(), Nt, MPI_DOUBLE,
                 nullptr, 0, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    }
    else if (task == 2) {
      // calculate linear dipole by reducing field strength by scale_factor and
      // multiplying response by same scale_factor
      double scale_factor = 1e3;
      for (auto& f: field_atomic) {
        f /= scale_factor;
      }
      argon.reset_to_ground_state();
      argon.calculate_dipole(field_atomic, field_atomic_dt, atomic_dt, linear_dipole);
      for (auto& d: linear_dipole) {
        d *= scale_factor;
      }
      MPI_Gather(linear_dipole.data(), Nt, MPI_DOUBLE,
                 nullptr, 0, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    }
    else {
      break;
    }
  }
}
