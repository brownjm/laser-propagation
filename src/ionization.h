#ifndef IONIZATION_H_
#define IONIZATION_H_

#include <vector>
#include <memory>
#include "interpolate.h"
#include "array2d.h"
#include "util.h"

class Radial;


namespace Ionization {


  class Rate {
  public:
    virtual double ionization_rate(double electric_field) = 0;
  };

  class TabulatedRate : public Rate {
  public:
    TabulatedRate(const std::string& filename, double scaling);

    double ionization_rate(double electric_field) override;
    
  private:
    std::string filename;
    double scaling;
    std::unique_ptr<Interpolate> interp;
  };
  


  class IonizationModel {
  public:
    IonizationModel(double density_of_neutrals, double ionizing_fraction, double pressure,
                    Rate& rate, int Nradius, int Ntime);
    void calculate_electron_density(const Radial& electric_field, Array2D<double>& electron_density);


    double density_of_neutrals, ionizing_fraction;
    Rate& rate;

    // for performance, to only calculate the rate once per step
    Array2D<double> cached_rate;

  };


}

#endif // IONIZATION_H_
