#include "electron_density.h"
#include "io.h"

namespace Observers {

  ElectronDensity::ElectronDensity(const std::string& fn)
    :filename(fn) {}
  
  void ElectronDensity::notify(int current_step, double,
		      const SimulationData& data) {
    IO::write_binary(IO::enumerate_filename(filename, current_step),
		     data.electron_density.vec());
  }
  
  void ElectronDensity::finalize() {}
}
