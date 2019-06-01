#include "parameters.h"
#include "propagator.h"
#include "driver.h"
#include "gaussian.h"
#include "fromfile.h"
#include "constants.h"
#include "medium.h"
#include "linear.h"
#include "results.h"
#include "timer.h"
#include "io.h"
#include "kerr.h"
#include "ramankerr.h"
#include "plasma.h"
#include "nonlinear_absorption.h"
#include "generate_rate.h"
#include "tabulated_rate.h"
#include "argon_response_pool.h"
#include "argon_response_worker.h"
#include <mpi.h>

template <typename T> std::string type_name();

void initialize_laser_field(Propagator& prop, Parameters::Parameters& params);
void initialize_linear_medium(Propagator& prop, Parameters::Parameters& params);
void initialize_results(Driver& driver, Parameters::Parameters& params);

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int rank, processors;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processors);
  if (rank == 0) {
    if (processors < 2) {
      std::cerr << "Number of processors (" << processors << ") must be 2 or greater\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
      return EXIT_FAILURE;
    }
    if (argc < 2) {
      std::cerr << "Please provide an input filename\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
      return EXIT_FAILURE;
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    try {
      Parameters::Parameters p(argv[1]);
      
      // coordinates
      int Ntime = p.get<int>("time/N");
      double time_min = p.get<double>("time/time_min");
      double time_max = p.get<double>("time/time_max");
      double wavelength_min = p.get<double>("time/wavelength_min");
      double wavelength_max = p.get<double>("time/wavelength_max");
      int Nradius = p.get<int>("space/N");

      if (Nradius != processors){
        std::cerr << "Number of points in radial grid (" << Nradius << ") must be equal to number of processors (" << processors << ")\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return EXIT_FAILURE;
      }
      
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
        double n0 = p.get<double>("calculated/n0");
        prop.add_polarization(std::make_shared<Kerr>(n2, n0, pressure));
      }
      else if (p.section_exists("ramankerr")) {
        double pressure = p.get<double>("medium/pressure");
        double n2 = p.get<double>("ramankerr/n2");
        double n0 = p.get<double>("calculated/n0");
        double fraction = p.get<double>("ramankerr/fraction");
        double gamma = p.get<double>("ramankerr/gamma");
        double lambda = p.get<double>("ramankerr/lambda");
        prop.add_polarization(std::make_shared<RamanKerr>(n2, n0, fraction, gamma, lambda, pressure));
      }
      
      if (p.section_exists("ionization")) {
        std::string filename;
        if (p.key_exists("ionization/read")) {
          filename = p.get<std::string>("ionization/read");
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
          IO::write(filename, intensities, rates);
        }
        else {
          throw std::runtime_error("Section [ionization] must contain either 'filename' or 'generate'");
        }
        
        double density_of_neutrals = p.get<double>("ionization/density_of_neutrals");
        double pressure = p.get<double>("medium/pressure");
        double fraction = p.get<double>("ionization/ionizing_fraction");
        
        auto ioniz = std::make_shared<TabulatedRate>(filename, density_of_neutrals, pressure,
                                                     fraction);
        prop.add_ionization(ioniz);
        
        double collision_time = p.get<double>("ionization/collision_time");
        prop.add_current(std::make_shared<Plasma>(collision_time, pressure));
        
        double ionization_potential = p.get<double>("ionization/ionization_potential");
        prop.add_current(std::make_shared<NonlinearAbsorption>(ionization_potential,
                                                               density_of_neutrals,
                                                               pressure, fraction,
                                                               prop.ionization_rate));
      }
      if (p.section_exists("argon")) {
        int Nr = p.get<int>("argon/Nr");
        int Nl = p.get<int>("argon/Nl");
        int Nmask = p.get<int>("argon/Nmask");
        std::string filename = p.get<std::string>("argon/potential");
        int ionization_box_size = p.get<int>("argon/ionization_box_size");
        double atomic_dt = p.get<double>("argon/step_size");
        int Nradius = p.get<int>("space/N");
        int Ntime = p.get<int>("time/N");
        double density_of_neutrals = p.get<double>("argon/density_of_neutrals");
        double pressure = p.get<double>("medium/pressure");
        auto argon = std::make_shared<ArgonResponsePool>(Nr, Nl, Nmask, filename,
                                                         ionization_box_size,
                                                         Nradius, Ntime, atomic_dt,
                                                         density_of_neutrals*pressure);
        
        // calculate ionization and nonlinear polarization using Argon
        prop.add_ionization(argon);
        prop.add_polarization(argon);
        
        // add plasma effects
        double collision_time = p.get<double>("argon/collision_time");
        prop.add_current(std::make_shared<Plasma>(collision_time, pressure));
        
        // add absorption from ionization
        double ionization_potential = p.get<double>("argon/ionization_potential");
        prop.add_current(std::make_shared<NonlinearAbsorption>(ionization_potential,
                                                               density_of_neutrals,
                                                               pressure, 1,
                                                               prop.ionization_rate));
      }
      
      Driver driver(prop);
      initialize_results(driver, p);
      
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
    stop_workers();
  }
  else {
    Parameters::Parameters p(argv[1]);
    
    // coordinates
    int Ntime = p.get<int>("time/N");
    double time_min = p.get<double>("time/time_min");
    double time_max = p.get<double>("time/time_max");
    double field_dt = (time_max - time_min) / (Ntime - 1);
    const double au_time = 2.418884326509e-17;
    double field_atomic_dt = field_dt / au_time;
    int Nr = p.get<int>("argon/Nr");
    int Nl = p.get<int>("argon/Nl");
    int Nmask = p.get<int>("argon/Nmask");
    std::string filename = p.get<std::string>("argon/potential");
    int ionization_box_size = p.get<int>("argon/ionization_box_size");
    double atomic_dt = p.get<double>("argon/step_size");
    ArgonResponseWorker argon(Nr, Nl, Nmask, filename, ionization_box_size,
                              Ntime, atomic_dt, field_atomic_dt);

    argon.run();
  }
  MPI_Finalize();
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

  if (p.section_exists("kerr")) {
    double n2 = p.get<double>("kerr/n2");
    double k0 = 2*Constants::pi / wavelength;
    double Pcr = 3.77*Constants::pi*n0 / (2*std::pow(k0, 2) * n2);
    p.set("calculated/Pcr", Pcr);
  }
  if (p.section_exists("ramankerr")) {
    double n2 = p.get<double>("ramankerr/n2");
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
void conditionally_add(Driver& driver, Parameters::Parameters& p, const std::string& name, ResultType result_type) {
  if (p.key_exists(name)) {
    std::string filename = p.get<std::string>(name);
    driver.add_result(std::make_unique<T>(filename), result_type);
  }
}

void initialize_results(Driver& driver, Parameters::Parameters& p) {
  // Add the coordinates and other diagnostics that are only run once
  conditionally_add<Results::Time>(driver, p, "results/time", ResultType::Once);
  conditionally_add<Results::Radius>(driver, p, "results/radius", ResultType::Once);
  conditionally_add<Results::Omega>(driver, p, "results/omega", ResultType::Once);
  conditionally_add<Results::Kperp>(driver, p, "results/kperp", ResultType::Once);
  conditionally_add<Results::Wavelength>(driver, p, "results/wavelength",
                                           ResultType::Once);
  conditionally_add<Results::Hankel>(driver, p, "results/hankel", ResultType::Once);
  conditionally_add<Results::Coef>(driver, p, "results/coef", ResultType::Once);
  conditionally_add<Results::Index>(driver, p, "results/index", ResultType::Once);
  conditionally_add<Results::Kz>(driver, p, "results/kz", ResultType::Once);

  // Add the field and density results
  conditionally_add<Results::Distance>(driver, p, "results/distance",
                                         ResultType::Expensive);
  conditionally_add<Results::TemporalField>(driver, p, "results/temporal_field",
                                              ResultType::Expensive);
  conditionally_add<Results::SpectralField>(driver, p, "results/spectral_field",
                                              ResultType::Expensive);
  conditionally_add<Results::ElectronDensity>(driver, p, "results/electron_density",
                                                ResultType::Expensive);

  // Add cheaper results
  conditionally_add<Results::Energy>(driver, p, "results/energy", ResultType::Cheap);
  conditionally_add<Results::MaxIntensity>(driver, p, "results/max_intensity",
                                             ResultType::Cheap);
  conditionally_add<Results::MaxDensity>(driver, p, "results/max_density",
                                           ResultType::Cheap);
}
