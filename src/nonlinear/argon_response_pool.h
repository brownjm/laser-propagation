#ifndef ARGON_RESPONSE_POOL_H_
#define ARGON_RESPONSE_POOL_H_

#include "ionization.h"
#include "nonlinear_response.h"
#include "argon.h"
#include "../core/array2d.h"

void prepare_workers_for_electron_density();
void prepare_workers_for_response();
void stop_workers();

class ArgonResponsePool : public Ionization, public NonlinearResponse {
public:
  ArgonResponsePool(int Nr, int Nl, int Nmask,
                const std::string& filename_potentials,
                double ionization_box_size,
                int Nradius, int Nt, double atomic_dt,
                double density_of_neutrals);

  void calculate_electron_density(const Radial& electric_field,
                                  Array2D<double>& ionization_rate,
                                  Array2D<double>& electron_density) override;
  
  void calculate_response(const std::vector<double>& radius,
                          const std::vector<double>& time,
                          const Array2D<std::complex<double>>& electric_field,
                          const Array2D<double>& electron_density,
                          Array2D<std::complex<double>>& response) override;
  
private:
  Argon argon;
  int Nradius, Nt;
  double atomic_dt, density_of_neutrals;
  const double au_efield = 514220670700.0;
  const double au_dipole = 8.47835309e-30;
  const double au_time = 2.418884326509e-17;
  Array2D<double> field_atomic, dipole, probability_free, dipole_linear;
  std::vector<double> local_field_atomic, local_dipole, local_probability_free, local_dipole_linear;
};

#endif // ARGON_RESPONSE_POOL_H_

