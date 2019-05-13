#ifndef SPECTRAL_FILTER_H_
#define SPECTRAL_FILTER_H_

#include "observer.h"

namespace Observers {

  class SpectralFilter : public Observer {
  public:
    SpectralFilter(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write("spectral_filter.dat", data.field.spectral_filter);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // SPECTRAL_FILTER_H_
