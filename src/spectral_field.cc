#include "spectral_field.h"
#include "io.h"

#include "observer.h"

namespace Observers {

  SpectralField::SpectralField(const std::string& fn)
    :filename(fn) {}
  
  void SpectralField::notify(int current_step, double, const SimulationData& data) {
    IO::write_binary(IO::enumerate_filename(filename, current_step),
		     data.field.spectral.vec());
  }
  
  void SpectralField::finalize() {}
}
