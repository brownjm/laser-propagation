#ifndef ARRAY2D_H_
#define ARRAY2D_H_

#include <complex>
#include <vector>

// a 2D interface to a vector

template <class T>
class Array2D {
public:
  Array2D() {}
  
  Array2D(int numrows, int numcols)
    :numrows(numrows), numcols(numcols), values(numrows*numcols) {}

  int rows() const { return numrows; }
  int cols() const { return numcols; }

  void resize(int numrows, int numcols) {
    this->numrows = numrows;
    this->numcols = numcols;
    values.resize(numrows*numcols);
  }
  
  // access to data, col is the fast axis
  inline T operator() (int row, int col) const {
    return values[row*numcols + col];
  }
  inline T& operator() (int row, int col) {
    return values[row*numcols + col];
  }

  T* get_data_ptr() { return values.data(); }
  const std::vector<T>& vec() const { return values; }


  int numrows, numcols;
  std::vector<T> values;
};


#endif // ARRAY2D_H_
