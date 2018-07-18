#include <iostream>
#include <iomanip>
#include "driver.h"
#include "propagator.h"
#include "observer.h"
#include "util.h"



Driver::Driver(Propagator& prop)
  :current_step(0), current_distance(0), propagator(prop) {}

void Driver::add_observer(std::unique_ptr<Observers::Observer> obs) {
  observers.push_back(std::move(obs));
}

void Driver::run(double distance, int steps) {
  // gives the observers the initial data
  notify_observers();

  std::cout << "z [m]    Energy [J]   Imax [W/m^2]   Rhomax [1/m^3]\n";
  auto data = propagator.get_data();
  std::cout << std::fixed << std::setprecision(3);
  std::cout<<  current_distance << "    ";
  std::cout << std::scientific << Util::energy(data.field) << "    ";
  std::cout << Util::max_intensity(data.field) << "      ";
  std::cout << Util::max_density(data.electron_density) << "\n";

  // advance the simulation forward
  for (int i = 1; i <= steps; ++i) {
    double zi = i * distance / steps;
    propagator.nonlinear_step(current_distance, zi);
    current_step = i;
    notify_observers();

    auto data = propagator.get_data();
    std::cout << std::fixed << std::setprecision(3);
    std::cout<<  current_distance << "    ";
    std::cout << std::scientific << Util::energy(data.field) << "    ";
    std::cout << Util::max_intensity(data.field) << "      ";
    std::cout << Util::max_density(data.electron_density) << "\n";
  }

  finalize();
}

void Driver::notify_observers() {
  auto data = propagator.get_data();
  for (auto& obs: observers) {
    obs->notify(current_step, current_distance, data);
  }
}

void Driver::finalize() {
  for (auto& obs: observers) {
    obs->finalize();
  }
}

