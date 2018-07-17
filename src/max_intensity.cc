#include "max_intensity.h"
#include "util.h"
#include "io.h"

namespace Observers {

  MaxIntensity::MaxIntensity(const std::string& fn)
    :filename(fn) {
    IO::clear_contents(filename);
  }
  
  void MaxIntensity::notify(int, double, const SimulationData& data) {
    double I = Util::max_intensity(data.field);
    IO::write_append(filename, I);
  }
  
  void MaxIntensity::finalize() {}
}
