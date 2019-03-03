# pde1d

[![Build Status](https://travis-ci.org/jgoldfar/pde1d.svg?branch=fixup-build)](https://travis-ci.org/jgoldfar/pde1d)

1D Partial Differential Equation Solver for MATLAB and Octave

pde1d solves systems of partial differential equations (PDE) in a single
spatial variable and time.
The input is mostly compatible with the MATLAB function pdepe. 
Many pdepe examples will work with pde1d with only small changes. 

However, pde1d contains several enhancements which make it substantially more powerful than pdepe.
Specifically, pde1d allows any number of ordinary differential equations (ODE) to be coupled to the system of PDE.
One use of these ODE, for example, is to allow more complex boundary conditions at the two ends of the PDE domain.
Another benefit of pde1d relative to pdepe is improved performance, particularly when many mesh points are required for a converged solution.
Third, pde1d allows advanced users to specify the order of the approximation functions in the spatial domain. 


Several examples and basic documentation are included.
An excellent introduction to solving PDE with the pdepe function is Professor Howard's note,
[Partial Differential Equations in MATLAB 7.0](http://www.math.tamu.edu/~phoward/m442/pdemat.pdf).
His examples, modified for pde1d can be found in the examples directory.

## Build Instructions

pde1d is written in C++, is built using the Make build system, and relies 
on the following third-party libraries.

* The IDA library from [Sundials](http://computation.llnl.gov/projects/sundials)
is used for solution of the differential-algebraic equations. Note that in order for the library to be compatible with Eigen (and the rest of this package), the index type should be set to `int32_t` at compile time; that is, set `SUNDIALS_INDEX_TYPE=INT32_T` on the `cmake` command line. The KLU linear solver should also be enabled by setting `KLU_ENABLE=ON`. Look to the build instructions in the `.travis.yml` file in this repository to see how this is done on Linux and MacOS.

* The KLU library from [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) is required.

* The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ matrix class library,
is used throughout the code.

* There is also a small dependency on the [Boost](http://www.boost.org/) C++ string package.

### Building on Linux

You'll need the following `apt` packages to build Sundials + PDE1D on Debian Stretch (or Ubuntu Xenial)

```shell
apt-get -qq -y --no-install-recommends install \
    build-essential \
    curl \
    # For curl
    libcurl4-openssl-dev \
    cmake \
    make \
    ca-certificates \
    python \
    libsuitesparse-dev \
    octave \
    liboctave-dev \
    libeigen3-dev \
    libboost-dev
```

The Makefile is tested on Linux, so the paths set there for Eigen, Octave, etc. hould be up-to-date as of March 2019, but the paths to Eigen, Suitesparse, Octave, and Sundials can all be changed on the command line.

After building Sundials as described above, the object files and the final `mex` file can be built by running

```shell
make objects pde1d.mex
```

### Building on MacOS

The Sundials + PDE1D build is tested on MacOS with the following [Homebrew](https://brew.sh/) dependencies:

```shell
brew install octave boost eigen suitesparse
```

After building Sundials as described above, the object files and the final `mex` file can be built by running

```shell
make objects pde1d.mex
```

### Installing PDE1D in Octave

After generating the `mex` file, you'll need to place it in the directory you intend to use it from, or somewhere else on your Octave/MATLAB Path, along with the corresponding `m` file.

To make it available to all users on a Linux or MacOS system, put those files in

```
/usr/local/share/octave/site/m/
```
