#ifndef IONIZATION_H_
#define IONIZATION_H_

#include <vector>
#include <memory>
#include "interpolate.h"

namespace Ionization {

  class Rate {
  public:
    Rate(const std::string& filename);
    double operator()(double intensity);

  private:
    std::vector<double> intensities, rates;
    std::unique_ptr<Util::Interpolate> interp;
  };


}

#endif // IONIZATION_H_
