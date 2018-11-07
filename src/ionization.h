#ifndef IONIZATION_H_
#define IONIZATION_H_

#include <vector>
#include <memory>
#include "interpolate.h"
#include "array2d.h"
#include "util.h"

class Radial;


namespace Ionization {
  
  class Ionization {
  public:
    virtual void calculate_electron_density(const Radial& electric_field,
                                            Array2D<double>& electron_density) = 0;
  };

  class Tabulated : public Ionization {
  public:
    Tabulated(const std::string& filename, double density_of_neutrals, double ionizing_fraction, double pressure, double scaling);
    void calculate_electron_density(const Radial& electric_field,
                                    Array2D<double>& electron_density) override;

  private:
    std::string filename;
    double density_of_neutrals, ionizing_fraction, pressure, scaling;
    std::unique_ptr<Interpolate> interp;
  };


}

#endif // IONIZATION_H_
