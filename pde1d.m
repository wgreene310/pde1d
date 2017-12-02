% pde1d Numerical solution of systems of partial differential equations (PDE)
% with two independent variables: a spatial dimension, x and time, t. A system
% of ordinary differential equations (ODE) may optionally be coupled with the
% PDE system.
%
% Usage:
% solution = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts)
% solution = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts,options)
% [solution,odeSolution] = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts,...
%                           odeFunc, odeIcFunc,xOde)
% [solution,odeSolution] = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts,...
%                           odeFunc, odeIcFunc,xOde,options)
%
% Equations and Boundary Conditions:
% The PDE to be solved are expressed in the following form:
% 
% c*Du/Dt = x^(-m)*D(x^m*f)/Dx + s
% D-           indicates a derivative
% m-           type of coordinate system defining the spatial dimension. If
%              zero, the system is rectangular Cartesian, if one, cylindrical, and if two,
%              spherical.
% c, f, s-     coefficients returned from the user-defined function, pdeFunc. 
%              In general they may be functions of x, t, u, and Du/Dx. The
%              c coefficient is often referred to as the mass matrix, f the flux,
%              and s, the source.
%
% Boundary conditions must be defined at both the left and right ends of the 
% spatial domain (i.e. at meshPts(1) and meshPts(end)). At each end, the
% boundary conditions are expressed in the following form:
%
% p + q*f = 0
% f-    flux defined in the PDE above.
% p, q- coefficients returned from the user-defined function, bcFunc.
%       p may be a function of x, t, and u. q may be a function of x and t.
%       Entries in the q vector that are zero at t=0, are assumed to be zero
%       at all future times. Entries in the q vector that are non-zero at t=0,
%       are assumed to be non-zero at all future times.
%
% The pdeFunc, icFunc, and bcFunc arguments to the pde1d function are function
% handles to user-defined functions that have the following forms:
% [c, f, s] = pdeFunc(x, t, u, DuDx)
% The function must return c, f, and s each of which has dimensions N x 1 where
% N is the number of PDE in the system. If any entry in the c-coefficient
% is zero at t=0, it is assumed to be zero for all future times. If any
% entry in the c-coefficient is non-zero at t=0, it is assumed to be non-zero
% at all future times.  
%
% When ODE are included as part of the system, pdeFunc takes the following form:
% [c, f, s] = pdeFunc(x, t, u, DuDx, v, vDot)
% A vector of ODE variables, v, and their derivatives with respect to time, vDot,
% are passed as additional arguments.
%
% [pLeft, qLeft, pRight, qRight] = bcFunc(xLeft, uLeft, xRight, uRight, t)
% The function must return pLeft, qLeft, pRight, and qRight each of which
% has dimensions N x 1 where N is the number of PDE in the system.
%
% When ODE are included as part of the system, bcFunc takes the following form:
% [pLeft, qLeft, pRight, qRight] = bcFunc(xLeft, uLeft, xRight, uRight, t, ...
%                                         v, vDot)
% A vector of ODE variables, v, and their derivatives with respect to time, vDot,
% are passed as additional arguments.
%
% u0 = icFunc(x)
% The complete specification of a PDE system requires that initial conditions
% (solution at t=0) be defined by the user. The function must return the 
% initial condition, u0, at spatial location x. u0 is a vector which has
% dimensions N x 1 where N is the number of PDE in the system. Formally,
% the boundary conditions returned from bcFunc and the initial conditions
% returned from icFunc should agree but this is not strictly required by
% pde1d.
%
% When ODE are included as part of the system, two additional functions must
% be defined.
% f=odeFunc(t,v,vDot,x,u,DuDx,flux, dudt, du2dxdt)
% The system of ODE is defined as a vector function,
% f(t,v,vDot,x,u,DuDx,flux, dudt, du2dxdt)=0, where
% t - time
% v - vector of ODE values
% vDot - derivative of ODE variables with respect to time
% x - vector of spatial locations where the ODE couple with the PDE variables
% u - values of the PDE variables at the x locations
% DuDx - derivatives of the PDE variables with respect to x, evaluated at
%        the x locations
% flux - values of the flux defined by the PDE definition evaluated at the x locations
% dudt - derivatives of the PDE variables with respect to time, evaluated at
%        the x locations
% du2dxdt - sectond derivatives of the PDE variables with respect to 
%           x and time, evaluated at the x locations
%
% v0=odeIcFunc()
% A vector of initial values of the ODE variables, v0, must be returned.
%
% meshPts- this argument to pde1d specifies the location of the mesh points
%          in the spatial dimension. This is a vector of x-locations where
%          the entries are strictly increasing.
%          The accuracy of the solution strongly
%          depends on having a sufficiently fine mesh and generally, the
%          user should try at least two different mesh densities to assess
%          convergence.
%
% timePts- vector of time points where the user would like to have the solution
%          returned. 
%
% options- name=value pairs to change default values for the solver.
%          RelTol=1e-3, relative tolerance for converged solution
%          AbsTol=1e-6, absolute tolerance for converged solution
%          Vectorized=false, if set to true, pdeFunc is called with a vector
%                     of x-values and is expected to return values of c, f, and s
%                     for all of these x-values. Setting this option to true
%                     substantially improves performance.
%          MaxSteps=10000, maximum number of time steps allowed
%
% solution- pde1d returns the solution of the system of PDE in a Mt x Mx x N
%           dimensioned matrix where Mt is the number of time points in the
%           timePts argument, Mx is the number of mesh points in the mesh
%           points argument, and N is the number of PDE in the system.

% Copyright (C) 2016-2017 William H. Greene
