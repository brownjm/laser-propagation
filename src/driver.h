#ifndef DRIVER_H_
#define DRIVER_H_

#include <vector>
#include <memory>

class Propagator;
namespace Observers {
  class Observer;
}

class Driver {
public:
  Driver(Propagator& propagator);
  void add_observer(std::unique_ptr<Observers::Observer> obs);
  void run(double start_distance, double end_distance, int steps);

private:
  int current_step;
  double solver_distance, physical_distance;
  
  Propagator& propagator;

  std::vector<std::unique_ptr<Observers::Observer>> observers;
  void notify_observers();
  void finalize(); // signals the observers to perform a final computation

  void print_runtime_data();
};


#endif // DRIVER_H_
