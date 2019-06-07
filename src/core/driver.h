#ifndef DRIVER_H_
#define DRIVER_H_

#include <vector>
#include <memory>

class Propagator;
namespace Results {
  class Result;
}

enum class ResultType {Once, Cheap, Expensive};

class Driver {
public:
  Driver(Propagator& propagator);
  void add_result(std::unique_ptr<Results::Result> obs, ResultType obstype);
  void run(double start_distance, double end_distance, int steps_cheap, int steps_expensive);

private:
  int current_step;
  double current_distance;
  
  Propagator& propagator;

  std::vector<std::unique_ptr<Results::Result>> once, cheap, expensive;
  void notify_once();
  void notify_cheap();
  void notify_expensive();
  void finalize(); // signals the observers to perform a final computation

  void print_runtime_data();
};


#endif // DRIVER_H_
