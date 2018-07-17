#include "temporal_field.h"
#include "io.h"

namespace Observers {

  TemporalField::TemporalField(const std::string& fn)
    :filename(fn) {}
  
  void TemporalField::notify(int current_step, double, const SimulationData& data) {
    IO::write_binary(IO::enumerate_filename(filename, current_step),
		     data.field.temporal.vec());
  }
  
  void TemporalField::finalize() {}
}
