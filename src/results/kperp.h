#ifndef KPERP_H_
#define KPERP_H_

#include "result.h"

namespace Results {

  class Kperp : public Result {
  public:
    Kperp(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.field.kperp);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // KPERP_H_
