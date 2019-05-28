#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "io.h"
#include "driver.h"
#include "propagator.h"
#include "result.h"
#include "util.h"
#include "timer.h"


Driver::Driver(Propagator& prop)
  :current_step(0), current_distance(0), propagator(prop) {}

void Driver::add_result(std::unique_ptr<Results::Result> obs, ResultType result_type) {
  if (result_type == ResultType::Once) {
    once.push_back(std::move(obs));
  } else if (result_type == ResultType::Cheap) {
    cheap.push_back(std::move(obs));
  } else if (result_type == ResultType::Expensive) {
    expensive.push_back(std::move(obs));
  }
}

void Driver::run(double start_distance, double stop_distance, int steps_cheap, int steps_expensive) {
  current_distance = start_distance;
  
  // runtime data header
  Timer timer;
  std::ostringstream ss;
  ss << "*** Runtime values ***\n";
  ss << "Started: " << timer.timestamp() << "\n";
  ss << "z [m]        Energy [J]   Imax [W/m^2]   Rhomax [1/m^3]\n";
  std::cout << ss.str();
  IO::write_append("log", ss.str());
  ss.str(std::string());
  ss.clear();

  print_runtime_data();

  // calculate distances
  double dz_cheap = (stop_distance - start_distance) / steps_cheap;
  double dz_expensive = (stop_distance - start_distance) / steps_expensive;
  std::vector<double> distances;
  for (int i = 0; i <= steps_cheap; ++i) {
    distances.push_back(i*dz_cheap + start_distance);
  }
  for (int i = 0; i <= steps_expensive; ++i) {
    distances.push_back(i*dz_expensive + start_distance);
  }

  // remove duplicates
  std::sort(std::begin(distances), std::end(distances));
  distances.erase(std::unique(std::begin(distances), std::end(distances)),
                  std::end(distances));
  
  // gives the observers the initial data
  notify_once();
  notify_cheap();
  notify_expensive();

  // advance the simulation forward
  for (std::size_t i = 1; i < distances.size(); ++i) {
    current_step = i;
    double z_next = distances[i];
    propagator.nonlinear_step(current_distance, z_next);

    // check if distance propagated is near to an interval where a
    // cheap or expensive diagnostic should be performed
    if (std::abs(std::remainder(current_distance-start_distance, dz_cheap)) < 1e-6) {
      notify_cheap();
    }
    if (std::abs(std::remainder(current_distance-start_distance, dz_expensive)) < 1e-6) {
      notify_expensive();
    }
        
    print_runtime_data();
  }

  finalize();

  // runtime footer
  ss << "Ended:   " << timer.timestamp() << "\n";
  ss << "Elapsed: " << timer.elapsed() << "\n";
  std::cout << ss.str();
  IO::write_append("log", ss.str());
}

void Driver::notify_once() {
  auto data = propagator.get_data();
  for (auto& obs: once) {
    obs->notify(current_step, current_distance, data);
  }
}

void Driver::notify_cheap() {
  auto data = propagator.get_data();
  for (auto& obs: cheap) {
    obs->notify(current_step, current_distance, data);
  }
}

void Driver::notify_expensive() {
  auto data = propagator.get_data();
  for (auto& obs: expensive) {
    obs->notify(current_step, current_distance, data);
  }
}

void Driver::finalize() {
  for (auto& obs: once) {
    obs->finalize();
  }
  for (auto& obs: cheap) {
    obs->finalize();
  }
  for (auto& obs: expensive) {
    obs->finalize();
  }
}

void Driver::print_runtime_data() {
  auto data = propagator.get_data();
  std::ostringstream ss;
  // ss << std::fixed << std::setprecision(3);
  ss << std::scientific << std::setprecision(3);
  ss << current_distance << "    ";
  ss << Util::energy(data.field) << "    ";
  ss << Util::max_intensity(data.field) << "      ";
  ss << Util::max_density(data.electron_density) << "\n";

  std::cout << ss.str();
  IO::write_append("log", ss.str());
}
