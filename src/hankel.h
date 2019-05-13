#ifndef HANKEL_H_
#define HANKEL_H_

#include "observer.h"

namespace Observers {

  class Hankel : public Observer {
  public:
    Hankel(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.field.dht.vec(), data.field.Nradius, data.field.Nradius);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // HANKEL_H_
