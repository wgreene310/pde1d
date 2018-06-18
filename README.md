# pde1d
1D Partial Differential Equation Solver for MATLAB and Octave

pde1d solves systems of partial differential equations (PDE) in a single
spatial variable and time. The input is mostly compatible with the MATLAB function pdepe. 
Many pdepe examples will work with pde1d with only small
changes. 

However, pde1d contains several enhancements which make it substantially
more powerful than pdepe. Specifically, pde1d allows any number of ordinary differential equations (ODE) to be coupled to the system of PDE. One use of these ODE, for example, is to allow more complex boundary conditions at the two ends of the PDE domain. Another benefit of pde1d relative to pdepe is improved performance, particularly when many mesh points are required for a converged solution. Third, pde1d allows advanced users to specify the order of the approximation functions in the spatial domain. 


Several examples and basic documentation are included.
An excellent introduction to solving PDE with the pdepe function is
Professor Howard's note,
[Partial Differential Equations in MATLAB 7.0](http://www.math.tamu.edu/~phoward/m442/pdemat.pdf). His examples, modified for pde1d, can be
found in the examples directory.

pde1d is written in C++, is built using the CMake build system, and relies 
on the following third-party libraries.
The IDA library from [Sundials] 
(http://computation.llnl.gov/projects/sundials)
is used for solution of the differential-algebraic equations.
The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ matrix class library,
is used throughout the code. There is also a small dependency on
the [Boost](http://www.boost.org/) C++ string package.
So all of these libraries are required to build the software
from scratch. The file, CMakeLists.txt, has been used to build versions
of pde1d.mex for Windows and Linux.