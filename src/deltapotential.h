#ifndef DELTAPOTENTIAL_H_
#define DELTAPOTENTIAL_H_

#include "nonlinear_response.h"
#include "atom.h"
//#include <omp.h>

class DeltaPotential : public NonlinearResponse {
public:
  DeltaPotential(double B, int Nt, double deltat, double density_of_neutrals)
    :B(B), Nt(Nt), dt(deltat/au_time), Na(density_of_neutrals) {
    // for (int i = 0; i < omp_num thres();) {
    //   atoms.push_back(atom);
    // }
    atom = std::make_unique<OneDAtom::OneDAtom>(B, Nt, dt);
    field.resize(Nt);
    current.resize(Nt);
  }

  void calculate_temporal_response(const Radial& electric_field, const Array2D<double>&,
                                   Radial& response) override {

    //#pragma omp parallel
    //{
    
    //#pragma omp for schedule(static)
    for (int i = 0; i < response.Nradius; ++i) {


      for (int j = 0; j < response.Ntime; ++j) {
        // convert to atomic field units
        field[j] = electric_field.rt(i, j).real() / au_field;
      }
      //atom->calculateNLCurrent(field.data(), current.data());
      atoms[thread_id]->calculateNLCurrent(field.data(), current.data());
      for (int j = 0; j < response.Ntime; ++j) {
        // convert to SI units
        response.rt(i, j) = current[j] * Na * au_dipole / au_time;
      }


      
    }
    //}
  }

private:
  const double au_time = 2.418884326509e-17;
  const double au_field = 5.142206707e11;
  const double au_dipole = 8.478353552e-30;
  double B;
  int Nt;
  double dt, Na;
  std::unique_ptr<OneDAtom::OneDAtom> atom;
  //std::vector<std::unique_ptr<OneDAtom::OneDAtom>> atoms;
  std::vector<double> field;
  std::vector<double> current;
};

#endif // DELTAPOTENTIAL_H_
