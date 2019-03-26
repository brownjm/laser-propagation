#ifndef COEF_H_
#define COEF_H_

#include "observer.h"

namespace Observers {

  class Coef : public Observer {
  public:
    Coef(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write("coef.dat", data.propagator.coef.vec(),
                  data.field.Nkperp, data.field.Nomega);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // COEF_H_
