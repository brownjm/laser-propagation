#ifndef TEMPORAL_FILTER_H_
#define TEMPORAL_FILTER_H_

#include "result.h"

namespace Results {

  class TemporalFilter : public Result {
  public:
    TemporalFilter(const std::string& filename)
      :filename(filename) {}
    void notify(int current_step, double, const SimulationData& data) override {
      if (current_step == 0) {
        IO::write("temporal_filter.dat", data.field.temporal_filter);
      }
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // TEMPORAL_FILTER_H_
