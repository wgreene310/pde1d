# pde1d
1D Partial Differential Equation Solver for MATLAB and Octave

pde1d solves systems of partial differential equations (PDE) in a single
spatial variable and time. 
The input is mostly compatible with the MATLAB function pdepe. 
Many pdepe examples will work with pde1d with only small changes. 

However, `pde1d` contains several enhancements which make it substantially
more powerful than pdepe. 
Specifically, `pde1d` allows any number of ordinary differential equations (ODE) to be coupled to the system of PDE. 
One use of these ODE, for example, is to allow more complex boundary conditions at the two ends of the PDE domain. 
Another benefit of pde1d relative to `pdepe` is improved performance, particularly when many mesh points are required 
for a converged solution. 
Third, pde1d allows advanced users to specify the order of the approximation functions in the spatial domain. 

Two capabilities of `pdepe` are not currently supported by `pde1d`.
`pde1d` does not allow complex coefficients. Also, when the PDE is defined in
a cylindrical or spherical coordinate system and the left end of the domain 
starts at zero, `pdepe` uses special approximation functions to account
for the singularity at this point; `pde1d` does not.

Several examples and basic documentation are included.
An excellent introduction to solving PDE with the pdepe function is
Professor Howard's note,
[Partial Differential Equations in MATLAB 7.0](http://www.math.tamu.edu/~phoward/m442/pdemat.pdf). 
His examples, modified for pde1d, can be
found in the examples directory.

## Build Instructions ##

`pde1d` is written in C++, is built using the CMake build system, and relies 
on the following third-party libraries.
The IDA library from [SUNDIALS](http://computation.llnl.gov/projects/sundials)
is used for solution of the differential-algebraic equations.
At least version 3 of Sundials is required.
The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) 
C++ matrix class library, is used throughout the code. 
By default, pde1d uses the Eigen sparse LU solver to solve the linear
algebraic system in SUNDIALS. However, if SUNDIALS has been built with
the optional KLU sparse solver from 
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html), 
`pde1d` can also be built to support this option. 
In this case the SuiteSparse libraries are required.
There is also a small dependency on
the [Boost](http://www.boost.org/) C++ string package.
At least version 3.0 of CMake is required.

The first step in building is to download and unpack the `pde1d` source
into a local directory. That local directory will be referred to as
`<path_to_pde1d_source>` below. A second directory, referred to as the
build directory, must then be created.

### Linux ###
When all of the required libraries are installed in "standard" locations,
`pde1d` can be built from the command line by invoking the following
commands from the build directory:

    cmake -DCMAKE_INSTALL_PREFIX=<somewhere_on_matlab_or_octave_path> \
      <path_to_pde1d_source>
	make

where `<somewhere_on_matlab_or_octave_path>` is a directory where `Matlab`
or `Octave` can find it.

By default, the `CMakeLists.txt` file assumes that `pde1d` is being built
for `Octave`. To build for `Matlab`, the `PSE` CMake variable must 
be set as shown below:

    cmake  -DPSE=Matlab \
     -DCMAKE_INSTALL_PREFIX=<somewhere_on_matlab_or_octave_path> \
       <path_to_pde1d_source>
	make

Only the latest Linux distributions have installable packages for the
required Sundials libraries. So it is often necessary to download the
source code from the Sundials site and follow their installation
procedure (also based on CMake) before proceeding with the `pde1d` build.
If Sundials is installed in a non-standard location, it can be selected
as follows in the `pde1d` build:

    cmake  -DPSE=Matlab \
     -DCMAKE_INSTALL_PREFIX=<somewhere_on_matlab_or_octave_path>  \
      -DSUNDIALS_ROOT=<path_to_sundials> <path_to_pde1d_source>
	make


### Windows ###
The easiest way to run CMake on Windows platforms is using
the graphical interface `cmake-gui` that is part of the Windows
CMake distribution. The source and build directories along with the necessary
CMake configuration variables can be set in the GUI. Because Windows
doesn't have standard locations for the installation of third-party
development libraries, the following CMake variables typically must
be set: `EIGEN_ROOT`, `SUNDIALS_ROOT`, `CMAKE_INSTALL_PREFIX`.
CMake supports several native build systems on Windows including MS Visual
Studio. After CMake generates the Visual Studio solution, it can be opened
and built from Visual Studio using the "Build Solution" menu selection.