#ifndef NONLINEAR_RESPONSE
#define NONLINEAR_RESPONSE

#include "radial.h"
#include "array2d.h"


// This is a base class for all nonlinear medium responses
class NonlinearResponse {
public:
  
  // Overload this function with the actual medium response
  // and fill in the values of reponse.rt(i, j)
  virtual void calculate_temporal_response(const Radial& electric_field, const Array2D<double>& electron_density, Radial& response) = 0;

  // Overload this function if the spectral response needs to be
  // modified after the temporal_to_spectral() transform
  virtual void finalize_spectral_response(Radial&) {}
};

#endif // NONLINEAR_RESPONSE
