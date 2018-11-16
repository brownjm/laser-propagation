#include <iostream>
#include <iomanip>
#include <sstream>
#include "io.h"
#include "driver.h"
#include "propagator.h"
#include "observer.h"
#include "util.h"
#include "timer.h"


Driver::Driver(Propagator& prop)
  :current_step(0), solver_distance(0), physical_distance(0), propagator(prop) {}

void Driver::add_observer(std::unique_ptr<Observers::Observer> obs) {
  observers.push_back(std::move(obs));
}

void Driver::run(double start_distance, double stop_distance, int steps) {
  physical_distance = start_distance;

  // log the output to a file
  std::ostringstream ss;
  //ss << propagator.log_grid_info();
  
  // time execution
  Timer timer;
  ss << "*** Runtime values ***\n";
  ss << "Started: " << timer.timestamp() << "\n";

  // gives the observers the initial data
  notify_observers();

  auto data = propagator.get_data();
  
  ss << "z [m]    Energy [J]   Imax [W/m^2]   Rhomax [1/m^3]\n";
  
  ss << std::fixed << std::setprecision(3);
  ss << physical_distance << "    ";
  ss << std::scientific << Util::energy(data.field) << "    ";
  ss << Util::max_intensity(data.field) << "      ";
  ss << Util::max_density(data.electron_density) << "\n";

  std::cout << ss.str();
  IO::write_append("log", ss.str());
  ss.str(std::string());
  ss.clear();

  // advance the simulation forward
  double dt = (stop_distance - start_distance) / steps;
  for (int i = 1; i <= steps; ++i) {
    double zi = i * dt;
    propagator.nonlinear_step(solver_distance, zi);
    current_step = i;
    physical_distance = solver_distance + start_distance;
    notify_observers();

    auto data = propagator.get_data();
    ss << std::fixed << std::setprecision(3);
    ss << physical_distance << "    ";
    ss << std::scientific << Util::energy(data.field) << "    ";
    ss << Util::max_intensity(data.field) << "      ";
    ss << Util::max_density(data.electron_density) << "\n";


    std::cout << ss.str();
    IO::write_append("log", ss.str());
    ss.str(std::string());
    ss.clear();
  }

  finalize();

  ss << "Ended:   " << timer.timestamp() << "\n";
  ss << "Elapsed: " << timer.elapsed() << "\n";
  std::cout << ss.str();
  IO::write_append("log", ss.str());
}

void Driver::notify_observers() {
  auto data = propagator.get_data();
  for (auto& obs: observers) {
    obs->notify(current_step, physical_distance, data);
  }
}

void Driver::finalize() {
  for (auto& obs: observers) {
    obs->finalize();
  }
}

