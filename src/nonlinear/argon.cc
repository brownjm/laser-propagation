#include "argon.h"
#include "../util/io.h"
#include <iostream>
#include <fstream>


Argon::Argon(int Nr, int Nl, int Nmask, const std::string& filename_potentials,
             double ionization_box_size)
  :Nr(Nr), dr(0.25), Nl(Nl), Nmask(Nmask),
   ionization_box_size(ionization_box_size),
   radius(Nr), gnd(Nl,Nr), psi(Nl,Nr), aux(Nl,Nr),
   V(Nl,Nr),
   alpha(Nr), beta(Nr), cl(Nl), mask(Nr), diag(Nr),
   sin(Nl,Nr), cos(Nl,Nr),
   loss_lmax(0) {

  // initialize coordinates
  for (int i = 0; i < Nr; ++i) {
    radius[i] = (i + 0.5) * dr;
  }

  for (int i = 0; i < Nr; ++i) {
    int j = i + 1;
    alpha[i] = -1.0/(2.0*dr*dr) * j*j / (j*j - 0.25);
    beta[i] = 1.0/(dr*dr)*(j*j - j + 0.5) / (j*j - j + 0.25);
  }
  
  // initialize V
  std::vector<int> momentum(Nl);
  for (int l = 0; l < Nl; ++l) {
    cl[l] = (l+1) / std::sqrt((2*l+1) * (2*l+3));
    momentum[l] = l;
  }

  // potential - long range and centrifugal
  for (int l = 0; l < Nl; ++l) {
    for (int i = 0; i < Nr; ++i) {
      V(l,i) = -1.0/radius[i] + l*(l+1) / (2*radius[i]*radius[i]);
    }
  }
  std::vector<double> potl0, potl1, potl2;
  read_potential(filename_potentials, potl0, potl1, potl2);
  for (std::size_t i = 0; i < potl0.size(); ++i) {
    V(0, i) += potl0[i];
    V(1, i) += potl1[i];
    V(2, i) += potl2[i];
  }
  for (int l = 3; l < Nl; ++l) {
    for (std::size_t i = 0; i < potl2.size(); ++i) {
      V(l, i) += potl2[i];
    }
  }

  // mask for r coordinates
  for (int i = 0; i < Nr-1; ++i) {
    if (i < (Nr - Nmask)) {
      mask[i] = 1.0;
    }
    else {
      mask[i] = std::pow(std::cos((i+Nr-Nmask)*M_PI / (2*Nmask)), 0.125);
    }
  }
  mask[Nr-1] = 0.0;
}

void Argon::find_ground_state(int steps, double dt, double loss) {
  // fill l=1 with hydrogen gnd state
  for (int i = 0; i < Nr; ++i) {
    psi(1, i) = 2*radius[i]*std::exp(-radius[i]);
  }

  // propagate in imaginary time
  for (int n = 1; n < steps; ++n) {
    stepH0(dt, loss);
    normalize();
  }
  // save current wavefunction to ground state
  gnd.values = psi.values;
}

void Argon::reset_to_ground_state() {
  psi.values = gnd.values;
}

void Argon::apply_Leven(const Array2D<complex>& in, Array2D<complex>& out) {
  complex imagi(0, 1);
  for (int l = 0; l < Nl-1; l+=2) {
    for (int i = 0; i < Nr; ++i) {
      out(l,i) = cos(l,i)*in(l,i) - imagi*sin(l,i)*in(l+1,i);
      out(l+1,i) = -imagi*sin(l,i)*in(l,i) + cos(l,i)*in(l+1,i);
    }
  }
  for (int i = 0; i < Nr; ++i) {
    out(Nl-1,i) = in(Nl-1,i);
  }
}

void Argon::apply_Lodd(const Array2D<complex>& in, Array2D<complex>& out) {
  complex imagi(0, 1);
  for (int l = 1; l < Nl-1; l+=2) {
    for (int i = 0; i < Nr; ++i) {
      out(l,i) = cos(l,i)*in(l,i) - imagi*sin(l,i)*in(l+1,i);
      out(l+1,i) = -imagi*sin(l,i)*in(l,i) + cos(l,i)*in(l+1,i);
    }
  }
  for (int i = 0; i < Nr; ++i) {
    out(0,i) = in(0,i);
  }
}

void Argon::step(double F, double dt) {
  // compute L matrix entries
  for (int l = 0; l < Nl; ++l) {
    for (int i = 0; i < Nr; ++i) {
      double arg = 0.5*dt*F*cl[l]*radius[i];
      cos(l, i) = std::cos(arg);
      sin(l, i) = std::sin(arg);
    }
  }

  apply_Leven(psi, aux);
  apply_Lodd(aux, psi);
  
  stepH0(dt, 0);
  
  apply_Lodd(psi, aux);
  apply_Leven(aux, psi);

  // zero out wavefunction that is in lmax
  for (int i = 0; i < Nr; ++i) {
    loss_lmax += dr*std::norm(psi(Nl-1,i));
    psi(Nl-1,i) = 0;
  }

  // apply mask in r
  for (int l = 0; l < Nl; ++l) {
    for (int i = 0; i < Nr; ++i) {
      psi(l, i) *= mask[i];
    }
  }
}

void Argon::calculate_dipole(const std::vector<double>& field, double field_dt, double atomic_dt, std::vector<double>& dipoles) {
  int steps = int(field_dt / atomic_dt + 1);
  atomic_dt = field_dt / steps;
  dipoles[0] = dipole();
  for (std::size_t n = 1; n < field.size(); ++n) {
    double dF = (field[n] - field[n-1]) / steps;
    double F = field[n-1];

    // make atomic steps
    for (int i = 0; i < steps; ++i) {
      F += dF;
      step(F, atomic_dt);
    }

    dipoles[n] = dipole();
  }
}

void Argon::calculate_ionization(const std::vector<double>& field, double field_dt, double atomic_dt, std::vector<double>& ionization) {
  int steps = int(field_dt / atomic_dt + 1);
  atomic_dt = field_dt / steps;
  ionization[0] = ionized();
  for (std::size_t n = 1; n < field.size(); ++n) {
    double dF = (field[n] - field[n-1]) / steps;
    double F = field[n-1];

    // make atomic steps
    for (int i = 0; i < steps; ++i) {
      F += dF;
      step(F, atomic_dt);
    }
    ionization[n] = ionized();
  }
}

void Argon::calculate_dipole_ionization(const std::vector<double>& field,
                                        double field_dt, double atomic_dt,
                                        std::vector<double>& dipoles,
                                        std::vector<double>& ionization) {
  int steps = int(field_dt / atomic_dt + 1);
  atomic_dt = field_dt / steps;
  dipoles[0] = dipole();
  ionization[0] = ionized();
  for (std::size_t n = 1; n < field.size(); ++n) {
    double dF = (field[n] - field[n-1]) / steps;
    double F = field[n-1];
    
    // make atomic steps
    for (int i = 0; i < steps; ++i) {
      F += dF;
      step(F, atomic_dt);
    }
    
    dipoles[n] = dipole();
    ionization[n] = ionized();
  }
}

void Argon::normalize() {
  // calculate the norm of the wavefunction: sum_l integral |psi_l|^2 dr
  double sum = 0;
  for (int l = 0; l < Nl; ++l) {
    for (int i = 0; i < Nr; ++i) {
      sum += std::norm(psi(l,i));
    }
  }
  double norm = sum * dr;
  double sqrt_norm = std::sqrt(norm);
  // normalize the wavefunction
  for (int l = 0; l < Nl; ++l) {
    for (int i = 0; i < Nr; ++i) {
      psi(l,i) /= sqrt_norm;
    }
  }
}

double Argon::ionized() {
  double sum = 0;
  for (int l = 0; l < Nl; ++l) {
    for (int i = 0; i < ionization_box_size; ++i) {
      sum += std::norm(psi(l,i));
    }
  }
  double norm = sum * dr;
  return 1 - norm;
}

double Argon::energy() {
  complex sum = 0;
  for (int l = 0; l < Nl; ++l) {
    sum += std::conj(psi(l,0)) * ((beta[0] + V(l,0))*psi(l,0) + alpha[0]*psi(l,1));
    for (int i = 1; i < Nr-1; ++i) {
      sum += std::conj(psi(l,i)) * (alpha[i-1]*psi(l,i-1) + (beta[i] + V(l,i))*psi(l,i) + alpha[i]*psi(l,i+1));
    }
    sum += std::conj(psi(l,Nr-1)) * (alpha[Nr-2]*psi(l,Nr-2) + (beta[Nr-1] + V(l,Nr-1))*psi(l,Nr-1));
  }
  return std::real(sum*dr);
}

double Argon::dipole() {
  double sum = 0.0;
  // l = 0 special cas
  for (int i = 0; i < Nr; ++i) {
    sum += std::real(std::conj(psi(0,i)) * radius[i] * cl[0]*psi(1,i));
  }

  // l > 0
  for (int l = 1; l < Nl-1; ++l) {
    for (int i = 0; i < Nr; ++i) {
      sum += std::real(std::conj(psi(l,i)) * radius[i] * (cl[l]*psi(l+1,i) + cl[l-1]*psi(l-1, i)));
    }
  }
  return -sum*dr;
}

double Argon::accel(double F) {
  double sum = 0.0;

  // the dipole acceleration is a(t) = -2*Re(<2|1>)
  // where |1> = (H0 z - z H0) |psi> => one = a + b
  // a = (T z - z T) |psi>
  // b = (V z - z V) |psi>
  // note: V includes short-range, long-range, and centrifugal
  // |2> = (H0 + F z) |psi> => two = c + d


  // l = 0
  complex a = dr * (cl[0]*(alpha[0]*psi(1,1)));
  complex b = -radius[0] * (-cl[0]*(V(0,0)-V(1,0))*psi(1,0));
  complex one = a + b;
  complex c = (beta[0] + V(0,0))*psi(0,0) + alpha[0]*psi(0,1);
  complex d = F * radius[0] * (cl[0]*psi(1,0));
  complex two = c + d;
  sum += 2*std::real(std::conj(two)*one);
  for (int i = 1; i < Nr-1; ++i) {
    a = dr * (cl[0]*(alpha[i]*psi(1,i+1) - alpha[i-1]*psi(1,i-1)));
    b = -radius[i] * (-cl[0]*(V(0,i)-V(1,i))*psi(1,i));
    one = a + b;
    c = alpha[i-1]*psi(0,i-1) + (beta[i] + V(0,i))*psi(0,i) + alpha[i]*psi(0,i+1);
    d = F * radius[i] * (cl[0]*psi(1,i));
    two = c + d;
    sum += 2*std::real(std::conj(two)*one);
  }
  a = dr * (cl[0]*(-alpha[Nr-2]*psi(1,Nr-2)));
  b = -radius[Nr-1] * (-cl[0]*(V(0,Nr-1)-V(1,Nr-1))*psi(1,Nr-1));
  one = a + b;
  c = alpha[Nr-2]*psi(0,Nr-2) + (beta[Nr-1] + V(0,Nr-1))*psi(0,Nr-1);
  d = F * radius[Nr-1] * (cl[0]*psi(1,Nr-1));
  two = c + d;
  sum += 2*std::real(std::conj(two)*one);

  // l > 1
  for (int l = 1; l < Nl-1; ++l) {
    a = dr * (cl[l-1]*(alpha[0]*psi(l-1,1)) + cl[l]*(alpha[0]*psi(l+1,1)));
    b = -radius[0] * (cl[l-1]*(V(l-1,0)-V(l,0))*psi(l-1,0) - cl[l]*(V(l,0)-V(l+1,0))*psi(l+1,0));
    one = a + b;
    c = (beta[0] + V(l,0))*psi(l,0) + alpha[0]*psi(l,1);
    d = F * radius[0] * (cl[l]*psi(l+1,0) + cl[l-1]*psi(l-1,0));
    two = c + d;
    sum += 2*std::real(std::conj(two)*one);
    for (int i = 1; i < Nr; ++i) {
      a = dr * (cl[l-1]*(alpha[i]*psi(l-1,i+1) - alpha[i-1]*psi(l-1,i-1)) + cl[l]*(alpha[i]*psi(l+1,i+1) - alpha[i-1]*psi(l+1,i-1)));
      b = -radius[i] * (cl[l-1]*(V(l-1,i)-V(l,i))*psi(l-1,i) - cl[l]*(V(l,i)-V(l+1,i))*psi(l+1,i));
      one = a + b;
      c = alpha[i-1]*psi(l,i-1) + (beta[i] + V(l,i))*psi(l,i) + alpha[i]*psi(l,i+1);
      d = F * radius[i] * (cl[l]*psi(l+1,i) + cl[l-1]*psi(l-1,i));
      two = c + d;
      sum += 2*std::real(std::conj(two)*one);
    }
    a = dr * (cl[l-1]*(-alpha[Nr-2]*psi(l-1,Nr-2)) + cl[l]*(-alpha[Nr-2]*psi(l+1,Nr-2)));
    b = -radius[Nr-1] * (cl[l-1]*(V(l-1,Nr-1)-V(l,Nr-1))*psi(l-1,Nr-1) - cl[l]*(V(l,Nr-1)-V(l+1,Nr-1))*psi(l+1,Nr-1));
    one = a + b;
    c = alpha[Nr-2]*psi(l,Nr-2) + (beta[Nr-1] + V(l,Nr-1))*psi(l,Nr-1);
    d = F * radius[Nr-1] * (cl[l]*psi(l+1,Nr-1) + cl[l-1]*psi(l-1,Nr-1));
    two = c + d;
    sum += 2*std::real(std::conj(two)*one);
  }
  return -sum*dr;
}

void Argon::save_wavefunction(const std::string& filename) {
  IO::write_binary(filename, psi.values);
}

void Argon::stepH0(double dt, double loss) {
  complex imagi(0, 1);
  complex cdt2 = 0.5*imagi*dt * (1.0 - imagi*loss);

  // apply [1 - idt/2 H0]
  for (int l = 0; l < Nl; ++l) {
    aux(l,0) = (1.0 - cdt2*(beta[0] + V(l,0)))*psi(l,0) - cdt2*alpha[0]*psi(l,1);
    for (int i = 1; i < Nr-1; ++i) {
      aux(l,i) = (1.0 - cdt2*(beta[i] + V(l,i)))*psi(l,i)
        -cdt2*(alpha[i-1]*psi(l,i-1) + alpha[i]*psi(l,i+1));
    }
    aux(l,Nr-1) = (1.0 - cdt2*(beta[Nr-1] + V(l,Nr-1)))*psi(l,Nr-1) - cdt2*alpha[Nr-2]*psi(l,Nr-2);
  }

  // solve [1 + idt/2 H0]^-1
  for (int l = 0; l < Nl; ++l) {
    // calculate diagonal [1 + idt/2 H0]
    for (int i = 0; i < Nr; ++i) {
      diag[i] = 1.0 + cdt2*(beta[i] + V(l,i));
    }

    // solve LHS
    for (int i = 1; i < Nr; ++i) {
      complex w = cdt2*alpha[i-1] / diag[i-1];
      diag[i] -= w * cdt2*alpha[i-1];
      aux(l,i) -= w * aux(l,i-1);
    }
    psi(l,Nr-1) = aux(l,Nr-1) / diag[Nr-1];
    for (int i = Nr-2; i >= 0; --i) {
      psi(l,i) = (aux(l,i) - cdt2*alpha[i]*psi(l,i+1)) / diag[i];
    }
  }
}


void read_potential(const std::string& filename,
                    std::vector<double>& potl0, std::vector<double>& potl1,
                    std::vector<double>& potl2) {
  std::ifstream input(filename);
  if (!input) throw std::runtime_error("Cannot open " + filename);
  int Nr;
  double dr;
  input >> Nr >> dr;
  double vl0, vl1, vl2;
  while (input >> vl0 >> vl1 >> vl2) {
    potl0.push_back(vl0);
    potl1.push_back(vl1);
    potl2.push_back(vl2);
  }
}


