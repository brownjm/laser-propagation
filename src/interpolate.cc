#include "interpolate.h"


#include <iomanip>
#include <stdexcept>

namespace Util {
  Interpolate::Interpolate(std::vector<double> xx, std::vector<double> yy)
    :xpts(xx), ypts(yy) {
    xmin = xpts.front();
    xmax = xpts.back();
    spline = gsl_spline_alloc(gsl_interp_linear, xpts.size());
    gsl_spline_init(spline, xpts.data(), ypts.data(), xpts.size());
    acc = gsl_interp_accel_alloc();
  }

  Interpolate::~Interpolate() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }

  double Interpolate::operator()(double x) {
    //return gsl_spline_eval(spline, x, acc);
    if ((x >= xmin) && (x <= xmax)) {
      return gsl_spline_eval(spline, x, acc);
    }
    else {
      std::ostringstream oss;
      oss << "Interpolation Error: ";
      oss << std::scientific << x << " is not within range (";
      oss << std::scientific << xmin << ", " << std::scientific << xmax << ")\n";
      throw std::runtime_error(oss.str());
    }
  }
}
