#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "observer.h"

namespace Observers {

  class Distance : public Observer {
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
