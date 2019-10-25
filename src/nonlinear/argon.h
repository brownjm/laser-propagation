#ifndef ARGON_H_
#define ARGON_H_

#include <complex>
#include <vector>
#include "../core/array2d.h"

typedef std::complex<double> complex;

class Argon {
public:
  Argon(int Nr, int Nl, int Nmask, const std::string& filename_potentials,
        double ionization_box_size);
  
  void find_ground_state(int steps, double dt, double loss);
  void reset_to_ground_state();

  void calculate_dipole(const std::vector<double>& field, double field_dt, double atomic_dt,
                        std::vector<double>& dipoles);

  void calculate_ionization(const std::vector<double>& field, double field_dt,
                            double atomic_dt, std::vector<double>& ionization);

  void calculate_dipole_ionization(const std::vector<double>& field,
                                   double field_dt, double atomic_dt,
                                   std::vector<double>& dipoles,
                                   std::vector<double>& ionization);
  
  double loss_in_lmax() {return loss_lmax;}
  double norm();
  void normalize();
  double energy();
  double dipole();
  double accel(double F);
  double ionized();
  void save_wavefunction(const std::string& filename);

private:
  int Nr;
  double dr;
  int Nl, Nmask;
  double ionization_box_size;
  std::vector<double> radius;
  Array2D<complex> gnd, psi, aux;
  Array2D<double> V;
  std::vector<double> alpha, beta, cl, mask;
  std::vector<complex> diag;
  Array2D<double> sin, cos;
  double loss_lmax;


  void stepH0(double dt, double loss=0);
  void apply_Leven(const Array2D<complex>& in, Array2D<complex>& out);
  void apply_Lodd(const Array2D<complex>& in, Array2D<complex>& out);
  void step(double F, double dt);  
};

void read_potential(const std::string& filename,
                    std::vector<double>& potl0, std::vector<double>& potl1,
                    std::vector<double>& potl2);

#endif // ARGON_H_
