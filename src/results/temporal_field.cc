#include "temporal_field.h"
#include "../util/io.h"

namespace Results {

  TemporalField::TemporalField(const std::string& fn)
    :number(0), filename(fn) {}
  
  void TemporalField::notify(int, double, const SimulationData& data) {
    IO::write_binary(IO::enumerate_filename(filename, number),
		     data.field.temporal.vec());
    number++;
  }
  
  void TemporalField::finalize() {}
}
