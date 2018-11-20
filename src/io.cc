#include <fstream>
#include <limits>
#include <iomanip>
#include "io.h"


void IO::read(const std::string& filename, std::vector<double>& x, std::vector<double>& y) {
  std::ifstream input(filename);
  if (!input) throw std::runtime_error("Cannot open input file: " + filename);
  double xi, yi;
  while (input >> xi >> yi) {
    x.push_back(xi);
    y.push_back(yi);
  }
}

void IO::write(const std::string& filename, const std::vector<double>& data) {
  std::ofstream output(filename);
  if (!output) throw std::runtime_error("Cannot open output file: " + filename);
  output.precision(std::numeric_limits<double>::max_digits10);
  output << std::scientific;
  for (auto& v: data) {
    output << v << "\n";
  }
}


void IO::write(const std::string& filename, const std::vector<std::complex<double>>& data,
	       int Nrows, int Ncols) {
  std::ofstream output(filename);
  if (!output) throw std::runtime_error("Cannot open output file: " + filename);
  output.precision(std::numeric_limits<double>::max_digits10);
  output << std::scientific;
  for (int r = 0; r < Nrows; ++r) {
    for (int c = 0; c < Ncols; ++c) {
      output << data[r*Ncols + c].real() << " " << data[r*Ncols + c].imag() << " ";
    }
    output << "\n";
  }
}

void IO::read_binary(const std::string& filename, std::vector<double>& data) {
  std::ifstream input(filename, std::ios::in | std::ios::binary);
  if (!input) throw std::runtime_error("Cannot open input file: " + filename);
  double value;
  while (input.read(reinterpret_cast<char*>(&value), sizeof(value))) {
    data.push_back(value);
  }
}

void IO::write_binary(const std::string& filename, const std::vector<double>& data) {
  std::ofstream output(filename, std::ios::out | std::ios::binary);
  if (!output) throw std::runtime_error("Cannot open output file: " + filename);
  output.write(reinterpret_cast<const char*>(data.data()),
	       data.size()*sizeof(double));
}

void IO::read_binary(const std::string& filename, std::vector<std::complex<double>>& data) {
  std::ifstream input(filename, std::ios::in | std::ios::binary);
  if (!input) throw std::runtime_error("Cannot open input file: " + filename);
  std::complex<double> value;
  while (input.read(reinterpret_cast<char*>(&value), sizeof(value))) {
    data.push_back(value);
  }
}

void IO::write_binary(const std::string& filename, const std::vector<std::complex<double>>& data) {
  std::ofstream output(filename, std::ios::out | std::ios::binary);
  if (!output) throw std::runtime_error("Cannot open output file: " + filename);
  output.write(reinterpret_cast<const char*>(data.data()),
	       data.size()*sizeof(std::complex<double>));
}


std::string IO::enumerate_filename(const std::string& name, int i) {
  std::ostringstream oss;
  oss << name << std::setfill('0') << std::setw(3) << i << ".dat";
  return std::move(oss.str());
}

void IO::clear_contents(const std::string& filename) {
  std::ofstream output(filename, std::ios::out | std::ios::trunc);
}

void IO::write_append(const std::string& filename, double value) {
  std::ofstream output(filename, std::ios::out | std::ios::app);
  if (!output) throw std::runtime_error("Cannot open output file: " + filename);
  output.precision(std::numeric_limits<double>::max_digits10);
  output << std::scientific;
  output << value << "\n";
}

void IO::write_append(const std::string& filename, const std::string& text) {
  std::ofstream output(filename, std::ios::out | std::ios::app);
  if (!output) throw std::runtime_error("Cannot open output file: " + filename);
  output << text;
}
