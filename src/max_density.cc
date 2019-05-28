#include "max_density.h"
#include "util.h"
#include "io.h"

namespace Results {

  MaxDensity::MaxDensity(const std::string& fn)
    :filename(fn) {
    IO::clear_contents(filename);
  }
  
  void MaxDensity::notify(int, double distance, const SimulationData& data) {
    double rhomax = Util::max_density(data.electron_density);
    IO::write_append(filename, distance, rhomax);
  }
  
  void MaxDensity::finalize() {}
}
