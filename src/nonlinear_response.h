#ifndef NONLINEAR_RESPONSE
#define NONLINEAR_RESPONSE

#include "radial.h"
#include "array2d.h"


// This is a base class for all nonlinear medium responses
class NonlinearResponse {
public:
  
  // Overload this function with the actual medium response and add
  // its values to response(i,j) where the indices (i, j) correspond
  // to (radius, time)
  // e.g. response(i,j) += func(electric_field(i,j).real());
  virtual void calculate(const std::vector<double>& radius,
                         const std::vector<double>& time,
                         const Array2D<std::complex<double>>& electric_field,
                         const Array2D<double>& electron_density,
                         Array2D<std::complex<double>>& response) = 0;
};

#endif // NONLINEAR_RESPONSE
