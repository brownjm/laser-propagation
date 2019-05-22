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
#include "ramankerr.h"
#include "plasma.h"
#include "nonlinear_absorption.h"

template <typename T> std::string type_name();

void initialize_laser_field(Propagator& prop, Parameters::Parameters& params);
void initialize_linear_medium(Propagator& prop, Parameters::Parameters& params);
void initialize_observers(Driver& driver, Parameters::Parameters& params);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please provide an input filename\n";
    return EXIT_FAILURE;
  }

  try {
    Parameters::Parameters p(argv[1]);

    // coordinates
    int Ntime = p.get<int>("time/N");
    double time_min = p.get<double>("time/time_min");
    double time_max = p.get<double>("time/time_max");
    double wavelength_min = p.get<double>("time/wavelength_min");
    double wavelength_max = p.get<double>("time/wavelength_max");
    int Nradius = p.get<int>("space/N");
    double radius_max = p.get<double>("space/radius_max");

    int Nkperp;
    std::string medium_type = p.get<std::string>("medium/type");
    if (medium_type == "capillary") {
      Nkperp = p.get<int>("capillary/modes");
    }
    else {
      Nkperp = Nradius;
    }

    double absolute_error = p.get<double>("ode/absolute_error");
    double relative_error = p.get<double>("ode/relative_error");
    double first_step = p.get<double>("ode/first_step");
    double starting_distance = p.get<double>("propagation/starting_distance");
    
    Propagator prop(Ntime, time_min, time_max, wavelength_min, wavelength_max,
		    Nradius, radius_max, Nkperp,
		    absolute_error, relative_error, first_step, starting_distance);




    double time_filter_min = p.get<double>("filtering/time_filter_min");
    double time_filter_max = p.get<double>("filtering/time_filter_max");
    double wavelength_filter_min = p.get<double>("filtering/wavelength_filter_min");
    double wavelength_filter_max = p.get<double>("filtering/wavelength_filter_max");
    prop.initialize_filters(time_filter_min, time_filter_max,
                            wavelength_filter_min, wavelength_filter_max);
    
    initialize_linear_medium(prop, p);
    initialize_laser_field(prop, p);
    if (p.section_exists("kerr")) {
      double pressure = p.get<double>("medium/pressure");
      double n2 = p.get<double>("kerr/n2");
      prop.add_polarization(std::make_unique<Kerr>(n2, pressure));
    }
    else if (p.section_exists("ramankerr")) {
      double pressure = p.get<double>("medium/pressure");
      double n2 = p.get<double>("ramankerr/n2");
      double fraction = p.get<double>("ramankerr/fraction");
      double gamma = p.get<double>("gamma");
      double lambda = p.get<double>("lambda");
      prop.add_polarization(std::make_unique<RamanKerr>(n2, fraction, gamma, lambda, pressure));
    }
    
    if (p.section_exists("ionization")) {
      std::string filename;
      if (p.key_exists("ionization/filename")) {
        filename = p.get<std::string>("ionization/filename");
      } else if (p.key_exists("ionization/generate")) {
        filename = p.get<std::string>("ionization/generate");
        std::string formula = p.get<std::string>("ionization/formula");
        double ionization_potential = p.get<double>("ionization/ionization_potential");
        double wavelength = p.get<double>("laser/wavelength");
        
        GenerateRate gen(ionization_potential, wavelength);
        std::function<double(double)> rate;
        if (formula == "adk") {
          rate = std::bind(&GenerateRate::adk, gen, std::placeholders::_1);
        }
        else if (formula == "mpi") {
          rate = std::bind(&GenerateRate::mpi, gen, std::placeholders::_1);
        }
        else if (formula == "tunnel") {
          rate = std::bind(&GenerateRate::tunnel, gen, std::placeholders::_1);
        }
        else if (formula == "yudin") {
          rate = std::bind(&GenerateRate::yudin, gen, std::placeholders::_1);
        }
        else if (formula == "ilkov") {
          rate = std::bind(&GenerateRate::ilkov, gen, std::placeholders::_1);
        }
        else {
          throw std::runtime_error("Unknown rate formula '" + formula + "'");
        }
        
        std::vector<double> intensities, rates;
        // add a rate to zero for zero intensity
        intensities.push_back(0);
        rates.push_back(0);
        // calc rate for 10^13 - 10^20 W/m^2
        for (int i = 13; i < 20; ++i) {
          for (int n = 100; n < 1000; ++n) {
            double I = 0.01*n*std::pow(10, i);
            intensities.push_back(I);
            rates.push_back(rate(I));
          }
        }

        // write out rate
        std::ofstream output(filename);
        for (std::size_t i = 0; i < intensities.size(); ++i) {
          output << intensities[i] << " " << rates[i] << "\n";
        }
      }
      else {
        throw std::runtime_error("Section [ionization] must contain either 'filename' or 'generate'");
      }
        
      double density_of_neutrals = p.get<double>("ionization/density_of_neutrals");
      double pressure = p.get<double>("medium/pressure");
      double fraction = p.get<double>("ionization/ionizing_fraction");

      auto ioniz = std::make_unique<Ionization>(filename, density_of_neutrals, pressure,
                                                fraction);
      prop.add_ionization(std::move(ioniz));
      prop.calculate_electron_density();


      double collision_time = p.get<double>("ionization/collision_time");
      prop.add_current(std::make_unique<Plasma>(collision_time, pressure));

      double ionization_potential = p.get<double>("ionization/ionization_potential");
      prop.add_current(std::make_unique<NonlinearAbsorption>(ionization_potential,
                                                             density_of_neutrals,
                                                             pressure, fraction,
                                                             prop.ionization_rate));
    }

    Driver driver(prop);
    initialize_observers(driver, p);
    
    double ending_distance = p.get<double>("propagation/ending_distance");
    int steps_cheap = p.get<int>("propagation/num_reports_cheap");
    int steps_expensive = p.get<int>("propagation/num_reports_expensive");

    std::stringstream ss;
    p.print(ss);
    ss << prop.log_grid_info();
    IO::clear_contents("log");
    IO::write_append("log", ss.str());
    std::cout << ss.str();

    driver.run(starting_distance, ending_distance, steps_cheap, steps_expensive);
  }
  catch (std::exception& err) {
    std::cerr << err.what() << "\n";
  }
}


void initialize_laser_field(Propagator& prop, Parameters::Parameters& p) {
  // check if this is a new simulation or a restart from a set of files
  std::string laser_type = p.get<std::string>("laser/type");
  if (laser_type == "gaussian") {
    double wavelength = p.get<double>("laser/wavelength");
    double waist = p.get<double>("laser/waist");
    double length = p.get<double>("laser/length");
    double energy = p.get<double>("laser/energy");
    double focus = p.get<double>("laser/focus");
    double chirp = p.get<double>("laser/chirp");
    double phase_deg = p.get<double>("laser/phase_deg");
    double phase = phase_deg / 180.0 * Constants::pi;
    double delay = p.get<double>("laser/delay");
    double gvd = p.get<double>("calculated/gvd");
    Field::Gaussian laser_field(wavelength, waist, focus, length, phase, delay, energy,
                                chirp, gvd);
    prop.initialize_field(laser_field);
  }
  else if (laser_type == "file") {
    std::string filename = p.get<std::string>("laser/filename");
    Field::FromFile laser_field(filename, prop.field.radius, prop.field.time);
    prop.initialize_field(laser_field);
  }
  else if (laser_type == "restart") {
    std::string spectral_file = p.get<std::string>("laser/filename");
    prop.restart_from(spectral_file);
  }
  else {
    throw std::runtime_error("Unsupported laser/type: " + laser_type);
  }

  // calculate Pin and Pcr
  double wavelength = p.get<double>("laser/wavelength");
  double waist = p.get<double>("laser/waist");
  double length = p.get<double>("laser/length");
  double energy = p.get<double>("laser/energy");
  double n0 = p.get<double>("calculated/n0");
  double tau = length*1.699/2;
  double I0 = energy / (std::pow(waist, 2) * std::pow(Constants::pi/2, 1.5) * tau);
  double Pin = Constants::pi * std::pow(waist, 2) * I0 / 2;
  p.set("calculated/Pin", Pin);

  if (p.section_exists("kerr") || p.section_exists("ramankerr")) {
    double n2 = p.get<double>("kerr/n2");
    double k0 = 2*Constants::pi / wavelength;
    double Pcr = 3.77*Constants::pi*n0 / (2*std::pow(k0, 2) * n2);
    p.set("calculated/Pcr", Pcr);
  }
}

void initialize_linear_medium(Propagator& prop, Parameters::Parameters& p) {
  std::string type = p.get<std::string>("medium/type");
  std::string index_name = p.get<std::string>("medium/index");
  double pressure = p.get<double>("medium/pressure");

  Medium::IndexFunction index = Medium::select_linear_index(index_name);
  double wavelength = p.get<double>("laser/wavelength");
  double omega0 = 2*Constants::pi*Constants::c / wavelength;
  auto n = index(omega0);
  p.set("calculated/n0", n.real());
  
  if (type == "capillary") {
    double radius_max = p.get<double>("space/radius_max");
    double cladding = p.get<double>("capillary/cladding");
    Linear::Capillary linear(index, radius_max, cladding, pressure);
    prop.initialize_linear(linear, omega0);

    double gvd = linear.gvd(0, omega0);
    p.set("calculated/gvd", gvd);
  }
  else if (type == "freespace") {
    Linear::FreeSpace linear(index);
    prop.initialize_linear(linear, omega0);
    
    double gvd = linear.gvd(0, omega0);
    p.set("calculated/gvd", gvd);
  }
  else if (type == "diffractionless") {
    Linear::DiffractionLess linear(index);
    prop.initialize_linear(linear, omega0);

    double gvd = linear.gvd(0, omega0);
    p.set("calculated/gvd", gvd);
  }
  
  else {
    throw std::runtime_error("Unsupported medium/type: " + type);
  }
}

template <class T>
void conditionally_add(Driver& driver, Parameters::Parameters& p, const std::string& name, ObserverType obstype) {
  if (p.key_exists(name)) {
    std::string filename = p.get<std::string>(name);
    driver.add_observer(std::make_unique<T>(filename), obstype);
  }
}

void initialize_observers(Driver& driver, Parameters::Parameters& p) {
  // Add the coordinates and other diagnostics that are only run once
  conditionally_add<Observers::Time>(driver, p, "output/time", ObserverType::Once);
  conditionally_add<Observers::Radius>(driver, p, "output/radius", ObserverType::Once);
  conditionally_add<Observers::Omega>(driver, p, "output/omega", ObserverType::Once);
  conditionally_add<Observers::Kperp>(driver, p, "output/kperp", ObserverType::Once);
  conditionally_add<Observers::Wavelength>(driver, p, "output/wavelength",
                                           ObserverType::Once);
  conditionally_add<Observers::Hankel>(driver, p, "output/hankel", ObserverType::Once);
  conditionally_add<Observers::Coef>(driver, p, "output/coef", ObserverType::Once);
  conditionally_add<Observers::Index>(driver, p, "output/index", ObserverType::Once);
  conditionally_add<Observers::Kz>(driver, p, "output/kz", ObserverType::Once);

  // Add the field and density observers
  conditionally_add<Observers::Distance>(driver, p, "output/distance",
                                         ObserverType::Expensive);
  conditionally_add<Observers::TemporalField>(driver, p, "output/temporal_field",
                                              ObserverType::Expensive);
  conditionally_add<Observers::SpectralField>(driver, p, "output/spectral_field",
                                              ObserverType::Expensive);
  conditionally_add<Observers::ElectronDensity>(driver, p, "output/electron_density",
                                                ObserverType::Expensive);

  // Add cheaper diagnostics
  conditionally_add<Observers::Energy>(driver, p, "output/energy", ObserverType::Cheap);
  conditionally_add<Observers::MaxIntensity>(driver, p, "output/max_intensity",
                                             ObserverType::Cheap);
  conditionally_add<Observers::MaxDensity>(driver, p, "output/max_density",
                                           ObserverType::Cheap);
}
