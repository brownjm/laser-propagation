#include <iostream>
#include "util/parameters.h"
#include "util/io.h"
#include "core/radial.h"
#include "nonlinear/argon_response.h"

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Please provide an input filename and an electric field filename\n";
    return EXIT_FAILURE;
  }

  try {
    Parameters::Parameters p(argv[1]);

    // read coordinates
    int Ntime = p.get<int>("time/N");
    double time_min = p.get<double>("time/time_min");
    double time_max = p.get<double>("time/time_max");
    double wavelength_min = p.get<double>("time/wavelength_min");
    double wavelength_max = p.get<double>("time/wavelength_max");

    // read in electric field file
    // convert E(r=0,t) into a Radial class
    std::vector<std::complex<double>> field_data;
    std::vector<double> onaxis(Ntime);
    IO::read_binary(argv[2], field_data);
    Radial electric_field(Ntime, time_min, time_max, wavelength_min, wavelength_max, 1, 1, 1);
    for (int j = 0; j < Ntime; ++j) {
      electric_field.rt(0, j) = field_data[j];
      onaxis[j] = std::real(field_data[j]);
    }
    std::vector<double> radius = electric_field.radius;
    std::vector<double> time = electric_field.time;

    IO::write("test_argon_onaxis_field.dat", time, onaxis);
    
    // set up argon atom
    int Nr = p.get<int>("argon/Nr");
    int Nl = p.get<int>("argon/Nl");
    int Nmask = p.get<int>("argon/Nmask");
    std::string filename = p.get<std::string>("argon/potential");
    int ionization_box_size = p.get<int>("argon/ionization_box_size");
    double atomic_dt = p.get<double>("argon/step_size");
    double density_of_neutrals = p.get<double>("argon/density_of_neutrals");
    double pressure = p.get<double>("medium/pressure");
    ArgonResponse argon(Nr, Nl, Nmask, filename, ionization_box_size, 1, Ntime, atomic_dt,
                        density_of_neutrals*pressure);

    // ionization
    Array2D<double> ionization_rate(1, Ntime);
    Array2D<double> electron_density(1, Ntime);
    Array2D<std::complex<double>> response(1, Ntime);

    // perform quantum calculations
    argon.calculate_electron_density(electric_field, ionization_rate, electron_density);
    argon.save_wavefunction("test_argon_wavefunction.dat");
    argon.calculate_response(radius, time, electric_field.temporal, electron_density, response);


    std::vector<double> nonlinear_dipole(Ntime);
    for (int j = 0; j < Ntime; ++j) {
      nonlinear_dipole[j] = std::real(response(0, j));
    }
    
    // output data
    IO::write("test_argon_ionization_rate.dat", time, ionization_rate.values);
    IO::write("test_argon_electron_density.dat", time, electron_density.values);
    IO::write("test_argon_nonlinear_dipole.dat", time, nonlinear_dipole);
  }
  catch (std::exception& err) {
    std::cerr << err.what() << "\n";
  }
}
