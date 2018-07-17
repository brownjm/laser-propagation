#include "ionization.h"
#include "io.h"

namespace Ionization {

  Rate::Rate(const std::string& filename) {
    IO::read(filename, intensities, rates);
    interp = std::make_unique<Util::Interpolate>(intensities, rates);
  }

  double Rate::operator()(double intensity) {
    return interp->operator()(intensity);
  }
  
}
