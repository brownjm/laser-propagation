#ifndef KPERP_H_
#define KPERP_H_

#include "observer.h"

namespace Observers {

  class Kperp : public Observer {
  public:
    Kperp(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write(filename, data.field.kperp);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // KPERP_H_
