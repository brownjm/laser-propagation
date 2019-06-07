#include "fromfile.h"
#include "../util/io.h"
#include <algorithm>

namespace Field {
  
  FromFile::FromFile(const std::string& filename, const std::vector<double>& radius,
                     const std::vector<double>& time)
    : filename(filename), radius(radius), time(time), Nradius(radius.size()),
      Ntime(time.size()) {
    IO::read_binary(filename, field);
}

  std::complex<double> FromFile::operator()(double, double r, double t) const {
    int i, j;
    auto iter_radius = std::find(std::cbegin(radius), std::cend(radius), r);
    if (iter_radius != std::cend(radius)) {
      i = std::distance(std::cbegin(radius), iter_radius);
    }
    else {
      throw std::runtime_error("Cannot find radius value: " + std::to_string(r));
    }

    auto iter_time = std::find(std::cbegin(time), std::cend(time), t);
    if (iter_time != std::cend(time)) {
      j = std::distance(std::cbegin(time), iter_time);
    }
    else {
      throw std::runtime_error("Cannot find time value: " + std::to_string(t));
    }
    
    return field[i*Ntime + j];
  }

}
