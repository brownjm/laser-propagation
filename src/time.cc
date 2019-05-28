#include "time.h"
#include "io.h"

namespace Results {

  Time::Time(const std::string& fn)
    :filename(fn) {}
  
  void Time::notify(int, double current_distance, const SimulationData&) {
    if (current_step == 0) {
      IO::write(fn_time, data.field.time);
    }
  }

  void Time::finalize() {}
}
