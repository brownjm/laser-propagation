#ifndef RESULT_H_
#define RESULT_H_


#include "simulation_data.h"

namespace Results {

class Result {
public:
  // all Result objects will receive these parameters at intervals
  virtual void notify(int current_step, double current_distance,
                      const SimulationData& data) = 0;

  // signals to the Result object to perform a final computation,
  // such as write their results to a file
  virtual void finalize() = 0;
};

}

#endif // RESULT_H_
