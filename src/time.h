#ifndef TIME_H_
#define TIME_H_

#include "observer.h"

namespace Observers {

  class Time : public Observer {
  public:
    Time(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write(filename, data.field.time);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // TIME_H_
