# pde1d
1D Partial Differential Equation Solver for MATLAB and Octave

pde1d solves systems of partial differential equations in a single
spatial variable and time. The input is mostly compatible with the MATLAB function pdepe. 
Many pdepe examples will work with pde1d with only small
changes. However, the underlying implementation is substantially 
different from pdepe. 


A basic Cmake file is provided that may be useful in building the
pde1d.mex file. Several examples
and basic documentation are included.

An excellent introduction to solving PDE with the pdepe function is
Professor Howard's note,
[Partial Differential Equations in MATLAB 7.0](http://www.math.tamu.edu/~phoward/m442/pdemat.pdf). His examples, modified for pde1d, can be
found in the examples directory.


pde1d relies on the IDA library from [Sundials] 
(http://computation.llnl.gov/projects/sundials)
for solution of the differential-algebraic equations.
The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ matrix class library,
is used throughout the code. There is also a small dependency on
the [Boost](http://www.boost.org/) C++ string package.
So all of these libraries are required to build the software
from scratch. The file, CMakeLists.txt, has been used to build versions
of pde1d.mex on the Windows platform.