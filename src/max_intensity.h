#ifndef MAX_INTENSITY_H_
#define MAX_INTENSITY_H_

#include "result.h"

namespace Results {

class MaxIntensity : public Result {
public:
  MaxIntensity(const std::string& filename);
  void notify(int current_step, double current_distance,
	      const SimulationData& data) override;
  void finalize() override;

private:
  std::string filename;
};

}

#endif // MAX_INTENSITY_H_
