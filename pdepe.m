function solution = pdepe(m, pdeFunc, icFunc, bcFunc, meshPts, timePts, options)
  % pdepe Numerical solution of systems of partial differential equations (PDE)
  % with two independent variables: a spatial dimension, x and time, t.
  %
  % This is a wrapper around pde1d provided for compatibility with MATLAB codes.
  %
  % Usage:
  % solution = pdepe(m,pdeFunc,icFunc,bcFunc,meshPts,timePts)
  % solution = pdepe(m,pdeFunc,icFunc,bcFunc,meshPts,timePts,options)
  %
  % For more usage details, see `pde1d`.
  % Copyright (C) 2019 Jonathan Goldfarb
  if ~exist('options', 'var')
    solution = pde1d(m, pdeFunc, icFunc, bcFunc, meshPts, timePts)
  else
    solution = pde1d(m, pdeFunc, icFunc, bcFunc, meshPts, timePts, options)
  end
end
