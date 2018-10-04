cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

add_executable(freespace.out
  src/freespace.cc
  src/driver.cc
  src/propagator.cc
  src/radial.cc
  src/io.cc
  src/gaussian.cc
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

target_compile_options(freespace.out PUBLIC -Ofast -Wall -Wextra -pedantic -DHAVE_INLINE)
target_link_libraries(freespace.out ${SOURCES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES})


add_executable(capillary.out
  src/capillary.cc
  src/driver.cc
  src/propagator.cc
  src/radial.cc
  src/io.cc
  src/gaussian.cc
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

target_compile_options(capillary.out PUBLIC -Ofast -Wall -Wextra -pedantic -DHAVE_INLINE)
target_link_libraries(capillary.out ${SOURCES} ${FFTW_LIBRARIES} ${GSL_LIBRARIES})
