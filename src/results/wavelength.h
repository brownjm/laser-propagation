#ifndef WAVELENGTH_H_
#define WAVELENGTH_H_

#include "result.h"

namespace Results {

  class Wavelength : public Result {
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
