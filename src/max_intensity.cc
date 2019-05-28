#include "max_intensity.h"
#include "util.h"
#include "io.h"

namespace Results {

  MaxIntensity::MaxIntensity(const std::string& fn)
    :filename(fn) {
    IO::clear_contents(filename);
  }
  
  void MaxIntensity::notify(int, double distance, const SimulationData& data) {
    double I = Util::max_intensity(data.field);
    IO::write_append(filename, distance, I);
  }
  
  void MaxIntensity::finalize() {}
}
