#ifndef FROM_FILE_H_
#define FROM_FILE_H_

#include <vector>
#include "field.h"

namespace Field {
  
  class FromFile : public Field {
  public:
    FromFile(const std::string& filename, const std::vector<double>& radius,
             const std::vector<double>& time);
    std::complex<double> operator()(double radius, double time) const override;

  private:
    const std::string& filename;
    const std::vector<double>& radius;
    const std::vector<double>& time;
    const int Nradius, Ntime;
    std::vector<std::complex<double>> field;
  };
}

#endif // FROM_FILE_H_
