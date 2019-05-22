# laser-propagation
This project is a C++ code for simulating the nonlinear propagation of ultrashort laser pulses in freespace or in a gas-filled capillary.

The model for propagation is a carrier-field resolved equation called the Unidirectional Pulse Propagation Equation (UPPE, M. Kolesik et al. Phys. Rev. Lett. 89, 283902 (2002)), where linear propagation and sources of nonlinearity are separately considered. The equation is solved using spectral methods and an adaptive-step RK45 solver.

The code is organized so that it is relatively straightforward to modify the linear propagation media (e.g. index of refraction), the nonlinearities, and the results that are calculated.

# Installation
The project uses CMake to compile the code and requires that you have installed the GNU Scientific Library (GSL), FFTW3, and Eigen3.

After downloading or cloning the repository, navigate to the root directory of the project and compile the project in a separate build directory using
```shell
mkdir build
cd build
cmake ..
make
```
This will produce the executable `main.out`

# Usage
## Input file
The type of simulation, inclusion of physical phenomena, and control of output files are controled with an input file containing the parameters of the simulation.

The input file syntax is similar to an ini file, in that pairs of key=value are grouped under sections (surrounded by square braces \[\]). Comments are denoted with a # and any text after them is ignored. Here is an example:
```shell
[propagation]
starting_distance = 0
ending_distance = 1 # final propagation distance of simulation
```

Since the input file is used to control the simulation there are some sections that are required and others that are optional. See the documentation for a complete description of all available sections.

## Running
 To run the simulation use
```shell
./main.out [filename]
```
During runtime, the simulation writes to the terminal and to a file called `log`. The file log contains a listing of the input parameters plus any calculated values:
```shell
*** Parameters ***
calculated/Pcr: 1.20036e+10
calculated/Pin: 4.6962e+10
filtering/time_filter_max: 150e-15
...
```

values related to the computational grids which hold the space-time and spectral values of the field:
```shell
*** Computational Grid ***
Supported Wavelengths: (1.50146e-07, 8.99377e-05)
Ntime    = 2048
Nomega   = 599
Nradius  = 200
Nkperp   = 200
ODE size = 119800 complex<double>
```

and finally, values of the field and execution times:
```shell
*** Runtime values ***
Started: 2019-05-10 17:12:27
z [m]        Energy [J]   Imax [W/m^2]   Rhomax [1/m^3]
0.000e+00    2.931e-04    5.564e+17      3.090e+22
2.000e-02    2.832e-04    1.023e+18      1.796e+23
4.000e-02    2.357e-04    1.096e+18      2.104e+23
6.000e-02    1.749e-04    1.097e+18      1.975e+23
8.000e-02    1.620e-04    9.268e+17      7.919e+22
1.000e-01    1.562e-04    9.700e+17      2.338e+23
Ended:   2019-05-10 17:13:15
Elapsed: 00:00:48
```

## Displaying the results
The directory `plotting` contains a set of Python3 scripts for displaying the results of a simulation run. These files require Python3, as well as the libraries numpy, scipy, and matplotlib to be installed.

Most scripts take two arguments: the input file path and a distance. For example, to plot the on-axis electric field at a distance of 0.25 m use:
```shell
python3 onaxis.py mysimulation/input_file 0.25
```

The provided scripts can display the following items: 
- space-time electric field
  - Ereal(r, t): temporal.py
  - Ereal(r=0,t): onaxis.py
  - Intensity(r, t): intensity.py
  - MaxIntensity(z): max_intensity.py
  - Energy(z): energy.py
- spectral representation of the electric field
  - angularly resolved spectrum |A(k,omega)|^2: spectral.py
  - angularly integrated spectrum |S(wavelength)|^2: spectrum.py
  - log |S(wavelength)|^2: log_spectrum.py
  - log |S(THz)|^2: thz_spectrum.py
- electron density
  - Rho(r, t): density.py
  - MaxRho(z): max_density.py


All of the plotting scripts use `load.py` to work with the data that was written during the simulation. If you wish to create your own plotting scripts, simply import the `load` module and use the class `Results` to load the data that you're interested in. For example,
```python
import load
results = load.Results(input_filename)
z, E = results.electric_field(distance)

(plotting code)
```
# Examples
There are two provided input files found in the directory `examples`, one for freespace propagation and one for propagation in a capillary.
