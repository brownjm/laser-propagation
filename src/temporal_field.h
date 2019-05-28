#ifndef TEMPORAL_FIELD_H_
#define TEMPORAL_FIELD_H_

#include "result.h"

namespace Results {

  class TemporalField : public Result {
  public:
    TemporalField(const std::string& filename);
    void notify(int current_step, double current_distance,
		const SimulationData& data) override;
    void finalize() override;
    
  private:
    int number;
    std::string filename;
  };
  
}

#endif // TEMPORAL_FIELD_H_
