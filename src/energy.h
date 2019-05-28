#ifndef ENERGY_H_
#define ENERGY_H_

#include <string>

#include "result.h"

namespace Results {

class Energy : public Result {
public:
  Energy(const std::string& filename);
  void notify(int current_step, double current_distance,
	      const SimulationData& data) override;
  void finalize() override;

private:
  std::string filename;
};

}

#endif // ENERGY_H_
