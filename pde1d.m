% pde1d Numerical solution of systems of partial differential equations (PDE)
% with two independent variables: a spatial dimension, x and time, t. 
% Usage:
% solution = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts)
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
% [pLeft, qLeft, pRight, qRight] = bcFunc(xLeft, uLeft, xRight, uRight, t)
% The function must return pLeft, qLeft, pRight, and qRight each of which
% has dimensions N x 1 where N is the number of PDE in the system.
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
% solution- pde1d returns the solution of the system of PDE in a Mt x Mx x N
%           dimensioned matrix where Mt is the number of time points in the
%           timePts argument, Mx is the number of mesh points in the mesh
%           points argument, and N is the number of PDE in the system.

%!demo
%!function heatCond
%!
%!x = linspace(0,1,10);
%!t = linspace(0,.05,40);
%!
%!pdeFunc = @(x,t,u,DuDx) heatpde(x,t,u,DuDx);
%!icFunc = @(x) heatic(x);
%!bcFunc = @(xl,ul,xr,ur,t) heatbc(xl,ul,xr,ur,t);
%!
%!m=0;
%!u = pde1d(m, pdeFunc,icFunc,bcFunc,x,t);
%!
%!figure; plot(t, u(:,end)); grid on;
%!xlabel('Time'); ylabel('Temperature');
%!title('Temperature at right end as a function of time');
%!figure; plot(x, u(end,:)); grid on;
%!xlabel('x'); ylabel('Temperature');
%!title('Temperature along the length at final time');
%!
%!end
%!
%!function [c,f,s] = heatpde(x,t,u,DuDx)
%!c = 1;
%!f = 10*DuDx;
%!s = 0;
%!end
%!
%!function u0 = heatic(x)
%!u0 = 0;
%!end
%!
%!function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
%!pl = ul-10;
%!ql = 0;
%!pr = 0;
%!qr = 1;
%!end

% Copyright (C) 2016 William H. Greene
