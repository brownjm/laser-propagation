#ifndef OMEGA_H_
#define OMEGA_H_

#include "observer.h"

namespace Observers {

  class Omega : public Observer {
  public:
    Omega(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.field.omega);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // OMEGA_H_
