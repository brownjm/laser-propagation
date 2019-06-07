#include "distance.h"
#include "../util/io.h"

namespace Results {

  Distance::Distance(const std::string& fn)
    :filename(fn) {
    IO::clear_contents(filename);
  }
  
  void Distance::notify(int, double current_distance, const SimulationData&) {
    IO::write_append(filename, current_distance);
  }

  void Distance::finalize() {}
}
