#ifndef IONIZATION_H_
#define IONIZATION_H_

#include <vector>
#include "array2d.h"

class Radial;

// base class for ionization models
class Ionization {
public:
  virtual void calculate_electron_density(const Radial& electric_field,
                                          Array2D<double>& ionization_rate,
                                          Array2D<double>& electron_density) = 0;
};

#endif // IONIZATION_H_
