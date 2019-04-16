#ifndef INDEX_H_
#define INDEX_H_

#include "observer.h"

namespace Observers {

  class Index : public Observer {
  public:
    Index(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.propagator.index, data.propagator.index.size(), 1);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // INDEX_H_
