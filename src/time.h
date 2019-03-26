#ifndef TIME_H_
#define TIME_H_

#include "observer.h"

namespace Observers {

  class Time : public Observer {
  public:
    Time(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write(filename, data.field.time);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // TIME_H_
