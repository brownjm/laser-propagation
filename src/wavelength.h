#ifndef WAVELENGTH_H_
#define WAVELENGTH_H_

#include "observer.h"

namespace Observers {

  class Wavelength : public Observer {
  public:
    Wavelength(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.field.wavelength);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // WAVELENGTH_H_
