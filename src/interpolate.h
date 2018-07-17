#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <vector>
#include <gsl/gsl_spline.h>

namespace Util {
  
  class Interpolate {
  public:
    Interpolate(std::vector<double> xpts, std::vector<double> ypts);
    ~Interpolate();

    double operator()(double x);

  private:
    std::vector<double> xpts, ypts;
    double xmin, xmax;
    gsl_spline* spline;
    gsl_interp_accel* acc;
  };

}

#endif // INTERPOLATE_H_
