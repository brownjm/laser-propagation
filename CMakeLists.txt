cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )


add_executable(main-${GIT_BRANCH}.out
  src/main.cc
  src/driver.cc
  src/propagator.cc
  src/radial.cc
  src/io.cc
  src/gaussian.cc
  src/fromfile.cc
  src/parameters.cc
  src/ionization.cc
  src/interpolate.cc
  src/util.cc
  src/medium.cc
  src/distance.cc
  src/energy.cc
  src/temporal_field.cc
  src/spectral_field.cc
  src/max_intensity.cc
  src/coordinates.cc
  src/max_density.cc
  src/electron_density.cc
  )

target_compile_options(main-${GIT_BRANCH}.out PUBLIC -Ofast -g -Wall -Wextra -pedantic -DHAVE_INLINE)
target_link_libraries(main-${GIT_BRANCH}.out ${SOURCES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES} Eigen3::Eigen)
