#ifndef FIELD_H_
#define FIELD_H_

#include <complex>

namespace Field {

  class Field {
  public:
    virtual std::complex<double> operator()(double z, double radius, double time) const = 0;
  };
  
}

#endif // FIELD_H_
