#ifndef ARGON_RESPONSE_WORKER_H_
#define ARGON_RESPONSE_WORKER_H_

#include <string>
#include "argon.h"

class ArgonResponseWorker {
public:
  ArgonResponseWorker(int Nr, int Nl, int Nmask, const std::string& filename_potentials,
                      double ionization_box, int Nt, double atomic_dt,
                      double field_atomic_dt);

  void run();

private:
  Argon argon;
  double Nt, atomic_dt, field_atomic_dt;
  std::vector<double> field_atomic, dipole, probability_free, linear_dipole;
};

#endif // ARGON_RESPONSE_WORKER_H_
