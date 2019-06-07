#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "result.h"

namespace Results {

  class Distance : public Result {
  public:
    Distance(const std::string& filename);
    void notify(int current_step, double current_distance,
		const SimulationData& data) override;
    void finalize() override;
    
  private:
    std::string filename;
  };
  
}

#endif // DISTANCE_H_
