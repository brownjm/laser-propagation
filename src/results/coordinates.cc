#include "coordinates.h"
#include "../util/io.h"

namespace Results {

  Coordinates::Coordinates(const std::string& filename_time,
			  const std::string& filename_radius,
			  const std::string& filename_omega,
			  const std::string& filename_kperp,
			  const std::string& filename_wavelength)
    :fn_time(filename_time), fn_radius(filename_radius), fn_omega(filename_omega),
     fn_kperp(filename_kperp), fn_wave(filename_wavelength) {}
  
  void Coordinates::notify(int current_step, double, const SimulationData& data) {
    if (current_step == 0) {
      IO::write(fn_time, data.field.time);
      IO::write(fn_radius,  data.field.radius);
      IO::write(fn_omega, data.field.omega);
      IO::write(fn_kperp, data.field.kperp);
      IO::write(fn_wave, data.field.wavelength);
    }
  }
  
  void Coordinates::finalize() {}
}
