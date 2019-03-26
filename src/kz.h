#ifndef KZ_H_
#define KZ_H_

#include "observer.h"

namespace Observers {

  class Kz : public Observer {
  public:
    Kz(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write("kz.dat", data.propagator.kz.vec(), data.field.Nkperp, data.field.Nomega);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // KZ_H_
