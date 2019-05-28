#include "energy.h"
#include "util.h"
#include "io.h"

namespace Results {

  Energy::Energy(const std::string& fn)
    :filename(fn) {
    IO::clear_contents(filename);
  }
  
  void Energy::notify(int, double distance, const SimulationData& data) {
    double energy = Util::energy(data.field);
    IO::write_append(filename, distance, energy);
  }
  
  void Energy::finalize() {}
}
