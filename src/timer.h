#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>

class Timer {
 public:
  Timer() {
    previous = std::chrono::system_clock::now();
  }

  std::string timestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
  }

  double elapsed_seconds() {
    auto now = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - previous).count();
    //auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - previous).count();

    previous = now;
    return elapsed;
  }

  std::string elapsed() {
    auto now = std::chrono::system_clock::now();
    auto diff = now - previous;
    auto h = std::chrono::duration_cast<std::chrono::hours>(diff);
    diff -= h;
    auto m = std::chrono::duration_cast<std::chrono::minutes>(diff);
    diff -= m;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(diff);
    diff -= s;
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(2) << h.count() << ":";
    ss << std::setfill('0') << std::setw(2) << m.count() << ":";
    ss << std::setfill('0') << std::setw(2) << s.count();
    return ss.str();
  }

 private:
  std::chrono::system_clock::time_point previous;
  std::stringstream ss;
};
