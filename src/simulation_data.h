#ifndef SIMULATION_DATA_H_
#define SIMULATION_DATA_H_

#include "radial.h"
#include "array2d.h"

struct SimulationData {
  const Radial& field;
  const Array2D<double>& electron_density;
};

#endif // SIMULATION_DATA_H_
