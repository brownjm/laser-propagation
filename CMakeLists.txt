cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# single threaded executable
add_executable(main.out
  src/main.cc
  src/core/driver.cc
  src/core/propagator.cc
  src/core/radial.cc
  src/util/io.cc
  src/util/parameters.cc
  src/util/interpolate.cc
  src/util/util.cc
  src/laser/gaussian.cc
  src/laser/fromfile.cc
  src/nonlinear/tabulated_rate.cc
  src/nonlinear/generate_rate.cc
  src/linear/medium.cc
  src/linear/linear.cc
  src/results/distance.cc
  src/results/energy.cc
  src/results/temporal_field.cc
  src/results/spectral_field.cc
  src/results/max_intensity.cc
  src/results/coordinates.cc
  src/results/max_density.cc
  src/results/electron_density.cc
  src/nonlinear/argon.cc
  src/nonlinear/argon_response.cc
  )

target_compile_options(main.out PUBLIC -Ofast -g -Wall -Wextra -pedantic -DHAVE_INLINE)
target_link_libraries(main.out ${SOURCES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES} Eigen3::Eigen)

# MPI executable
find_package(MPI REQUIRED)
add_executable(main-mpi.out
  src/main-mpi.cc
  src/core/driver.cc
  src/core/propagator.cc
  src/core/radial.cc
  src/util/io.cc
  src/util/parameters.cc
  src/util/interpolate.cc
  src/util/util.cc
  src/laser/gaussian.cc
  src/laser/fromfile.cc
  src/nonlinear/tabulated_rate.cc
  src/nonlinear/generate_rate.cc
  src/linear/medium.cc
  src/linear/linear.cc
  src/results/distance.cc
  src/results/energy.cc
  src/results/temporal_field.cc
  src/results/spectral_field.cc
  src/results/max_intensity.cc
  src/results/coordinates.cc
  src/results/max_density.cc
  src/results/electron_density.cc
  src/nonlinear/argon.cc
  src/nonlinear/argon_response_pool.cc
  src/nonlinear/argon_response_worker.cc
  )

target_compile_options(main-mpi.out PUBLIC -Ofast -g -Wall -Wextra -pedantic -DHAVE_INLINE)
target_link_libraries(main-mpi.out ${SOURCES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES}
  Eigen3::Eigen ${MPI_LIBRARIES} MPI::MPI_CXX)
