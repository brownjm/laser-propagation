#ifndef COEF_H_
#define COEF_H_

#include "result.h"

namespace Results {

  class Coef : public Result {
  public:
    Coef(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write("coef.dat", data.propagator.coef.vec(),
                data.field.Nkperp, data.field.Nomega);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // COEF_H_
