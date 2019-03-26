#ifndef RADIUS_H_
#define RADIUS_H_

#include "observer.h"

namespace Observers {

  class Radius : public Observer {
  public:
    Radius(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write(filename, data.field.radius);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // RADIUS_H_
