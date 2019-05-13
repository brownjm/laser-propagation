#ifndef SPECTRAL_FIELD_H_
#define SPECTRAL_FIELD_H_

#include "observer.h"

namespace Observers {

  class SpectralField : public Observer {
  public:
    SpectralField(const std::string& filename);
    void notify(int current_step, double current_distance,
		const SimulationData& data) override;
    void finalize() override;
    
  private:
    int number;
    std::string filename;
  };
  
}

#endif // SPECTRAL_FIELD_H_
