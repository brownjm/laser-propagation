#ifndef OMEGA_H_
#define OMEGA_H_

#include "observer.h"

namespace Observers {

  class Omega : public Observer {
  public:
    Omega(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write(filename, data.field.omega);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // OMEGA_H_
