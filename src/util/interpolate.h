#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <vector>
#include <iomanip>
#include <gsl/gsl_spline.h>

std::runtime_error generate_interpolation_error(double x, double xmin, double xmax);
  
class Interpolate {
public:
  Interpolate(std::vector<double> xpts, std::vector<double> ypts);
  ~Interpolate();
  
  
  double operator()(double x) {
    if ((x >= xmin) && (x <= xmax)) return gsl_spline_eval(spline, x, acc);
    else throw generate_interpolation_error(x, xmin, xmax);
  }
  
private:
  std::vector<double> xpts, ypts;
  double xmin, xmax;
  gsl_spline* spline;
  gsl_interp_accel* acc;
};



#endif // INTERPOLATE_H_
