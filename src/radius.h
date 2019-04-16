#ifndef RADIUS_H_
#define RADIUS_H_

#include "observer.h"

namespace Observers {

  class Radius : public Observer {
  public:
    Radius(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.field.radius);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // RADIUS_H_
