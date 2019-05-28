#include "electron_density.h"
#include "io.h"

namespace Results {

  ElectronDensity::ElectronDensity(const std::string& fn)
    :number(0), filename(fn) {}
  
  void ElectronDensity::notify(int, double, const SimulationData& data) {
    IO::write_binary(IO::enumerate_filename(filename, number),
		     data.electron_density.vec());
    number++;
  }
  
  void ElectronDensity::finalize() {}
}
