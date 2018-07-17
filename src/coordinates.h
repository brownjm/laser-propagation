#ifndef COORDINATES_H_
#define COORDINATES_H_

#include "observer.h"

namespace Observers {

  class Coordinates : public Observer {
  public:
    Coordinates(const std::string& filename_time,
		const std::string& filename_radius,
		const std::string& filename_omega,
		const std::string& filename_kperp,
		const std::string& filename_wavelength);
    
    void notify(int current_step, double current_distance,
		const SimulationData& data) override;
    void finalize() override;
    
  private:
    std::string fn_time, fn_radius, fn_omega, fn_kperp, fn_wave;
  };
  
}

#endif // COORDINATES_H_
