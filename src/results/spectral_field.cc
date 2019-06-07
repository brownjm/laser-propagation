#include "spectral_field.h"
#include "../util/io.h"

#include "result.h"

namespace Results {

  SpectralField::SpectralField(const std::string& fn)
    :number(0), filename(fn) {}
  
  void SpectralField::notify(int, double, const SimulationData& data) {
    IO::write_binary(IO::enumerate_filename(filename, number),
		     data.field.spectral.vec());
    number++;
  }
  
  void SpectralField::finalize() {}
}
