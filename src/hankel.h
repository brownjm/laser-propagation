#ifndef HANKEL_H_
#define HANKEL_H_

#include "observer.h"

namespace Observers {

  class Hankel : public Observer {
  public:
    Hankel(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write(filename, data.field.dht.vec(), data.field.Nradius, data.field.Nradius);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // HANKEL_H_
