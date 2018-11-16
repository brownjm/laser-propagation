#include "parameters.h"
#include "propagator.h"
#include "driver.h"
#include "gaussian.h"
#include "fromfile.h"
#include "constants.h"
#include "medium.h"
#include "linear.h"
#include "observers.h"
#include "timer.h"
#include "io.h"
#include "kerr.h"
#include "plasma.h"
#include "nonlinear_absorption.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please provide an input filename\n";
    return EXIT_FAILURE;
  }

  try {
    Parameters::Parameters p(argv[1]);

    // check which simulation is requested
    std::string sim_type = p.get<std::string>("simulation/type");

    
    // coordinates
    int Ntime = p.get<int>("time/N");
    double time_min = p.get<double>("time/time_min");
    double time_max = p.get<double>("time/time_max");
    double wave_min = p.get<double>("time/wave_min");
    double wave_max = p.get<double>("time/wave_max");
    int Nradius = p.get<int>("space/N");
    double R = p.get<double>("space/R");

    int Nkperp;
    if (sim_type == "capillary") {
      Nkperp = p.get<int>("capillary/modes");
    }
    else if (sim_type == "freespace") {
      Nkperp = Nradius;
    }
    else {
      throw std::runtime_error("Unsupported simulation: " + sim_type + "\nOptions are capillary or freespace");
    }

    double abserr = p.get<double>("ode/abserr");
    double relerr = p.get<double>("ode/relerr");
    double first_step = p.get<double>("ode/first_step");

    Propagator prop(Ntime, time_min, time_max, wave_min, wave_max,
		    Nradius, R, Nkperp,
		    abserr, relerr, first_step);


    // check if this is a new simulation or a restart from a set of files
    std::string restart = p.get<std::string>("simulation/starting_mode");
    double wavelength = p.get<double>("laser/wavelength");
    if (restart == "new") {
      // input field
      std::string laser_type = p.get<std::string>("laser/type");
      if (laser_type == "gaussian") {
	double length = p.get<double>("laser/length");
	double waist = p.get<double>("laser/waist");
	double focus = p.get<double>("laser/focus");
	double energy = p.get<double>("laser/energy");
	double phase_deg = p.get<double>("laser/phase_deg");
	double delay = p.get<double>("laser/delay");
	double phase = phase_deg / 180.0 * Constants::pi;
	Field::Gaussian laser_field(wavelength, waist, focus, length, phase, delay, energy);
	prop.initialize_field(laser_field);
      }
      else if (laser_type == "file") {
	std::string filename = p.get<std::string>("laser/filename");
	Field::FromFile laser_field(filename, prop.field.radius, prop.field.time);
	prop.initialize_field(laser_field);
      }
      else {
	throw std::runtime_error("Unsupported laser/type: " + laser_type);
      }
    }
    else if (restart == "resume") {
      std::string spectral_file = p.get<std::string>("simulation/spectral_file");
      prop.restart_from(spectral_file);
    }
    else {
      throw std::runtime_error("Unsupported simulation/starting_mode: " + restart);
    }

    
    std::string medium_name = p.get<std::string>("medium/medium");
    double pressure = p.get<double>("medium/pressure");
    double n2 = p.get<double>("medium/n2");
    double collision_time = p.get<double>("medium/collision_time");
    std::string rate_filename = p.get<std::string>("medium/rate");
    double fraction = p.get<double>("medium/fraction");
    double scaling = p.get<double>("medium/scaling");
    double density_of_neutrals = p.get<double>("medium/density_of_neutrals");
    double ionization_potential = p.get<double>("medium/ionization_potential");


    auto index = Medium::select_linear_index(medium_name);
    double omega0 = 2*Constants::pi*Constants::c / wavelength;
    if (sim_type == "capillary") {
      double radius = p.get<double>("capillary/radius");
      double cladding = p.get<double>("capillary/cladding");
      Capillary linear(index, radius, cladding, pressure);
      prop.initialize_linear(linear, omega0);
    }
    else if (sim_type == "freespace") {
      FreeSpace linear(index);
      prop.initialize_linear(linear, omega0);
    }

    Ionization::TabulatedRate rate(rate_filename, scaling);
    auto ionization_model = std::make_shared<Ionization::IonizationModel>(density_of_neutrals, fraction, pressure, rate, Ntime, Nradius);
    prop.add_ionization(ionization_model);

    prop.calculate_electron_density();
    
    prop.add_polarization(std::make_unique<Kerr>(n2, pressure));
    prop.add_current(std::make_unique<Plasma>(collision_time, pressure));
    prop.add_current(std::make_unique<NonlinearAbsorption>(ionization_potential, density_of_neutrals, pressure, fraction, ionization_model));


    double time_filter_min = p.get<double>("time/time_filter_min");
    double time_filter_max = p.get<double>("time/time_filter_max");
    double wave_filter_min = p.get<double>("time/wave_filter_min");
    double wave_filter_max = p.get<double>("time/wave_filter_max");
    prop.initialize_filters(time_filter_min, time_filter_max,
                            wave_filter_min, wave_filter_max);


    
    // output files
    std::string fn_temporal = p.get<std::string>("output/temporal_field");
    std::string fn_spectral = p.get<std::string>("output/spectral_field");
    std::string fn_density = p.get<std::string>("output/electron_density");
    std::string fn_dist = p.get<std::string>("output/distances");
    std::string fn_time = p.get<std::string>("output/time");
    std::string fn_radius = p.get<std::string>("output/radius");
    std::string fn_omega = p.get<std::string>("output/omega");
    std::string fn_kperp = p.get<std::string>("output/kperp");
    std::string fn_wave = p.get<std::string>("output/wavelength");
    std::string fn_energy = p.get<std::string>("output/energy");
    std::string fn_intensity = p.get<std::string>("output/max_intensity");
    std::string fn_maxdensity = p.get<std::string>("output/max_density");
    
    Driver driver(prop);
    driver.add_observer(std::make_unique<Observers::Distance>(fn_dist));
    driver.add_observer(std::make_unique<Observers::Coordinates>(fn_time, fn_radius,
    								 fn_omega, fn_kperp,
    								 fn_wave));
    driver.add_observer(std::make_unique<Observers::TemporalField>(fn_temporal));
    driver.add_observer(std::make_unique<Observers::SpectralField>(fn_spectral));
    driver.add_observer(std::make_unique<Observers::ElectronDensity>(fn_density));
    driver.add_observer(std::make_unique<Observers::Energy>(fn_energy));
    driver.add_observer(std::make_unique<Observers::MaxIntensity>(fn_intensity));
    driver.add_observer(std::make_unique<Observers::MaxDensity>(fn_maxdensity));

    
    double start_distance = p.get<double>("propagation/start_distance");
    double end_distance = p.get<double>("propagation/end_distance");
    int steps = p.get<int>("propagation/steps");

    std::stringstream ss;
    //p.print(ss);
    IO::clear_contents("log");
    IO::write_append("log", ss.str());
    std::cout << ss.str();
    driver.run(start_distance, end_distance, steps);
  }
  catch (std::exception& err) {
    std::cerr << err.what() << "\n";
  }
}
