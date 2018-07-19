// Copyright (C) 2016-2017 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <MatrixTypes.h>

class PDESolution;

class PDE1dDefn
{
public:
  PDE1dDefn();
  virtual int getNumPDE() const = 0;
  virtual int getCoordSystem() const { return 0; }
  virtual void evalIC(double x, RealVector &ic) = 0;
  virtual void evalODEIC(RealVector &ic) = 0;
  struct BC {
    RealVector pl, ql, pr, qr;
  };
  virtual void evalBC(double xl, const RealVector &ul,
    double xr, const RealVector &ur, double t, 
    const RealVector &v, const RealVector &vDot, BC &bc) = 0;
  struct PDE {
    RealVector c, f, s;
  };
  virtual void evalPDE(double x, double t, 
    const RealVector &u, const RealVector &DuDx, 
    const RealVector &v, const RealVector &vDot, PDE &pde) = 0;
  virtual const RealVector &getMesh() const = 0;
  virtual const RealVector &getTimeSpan() const = 0;
  virtual bool hasVectorPDEEval() const { return false;  }
  struct PDECoeff {
    RealMatrix c, f, s;
  };
  virtual void evalPDE(const RealVector &x, double t,
    const RealMatrix &u, const RealMatrix &DuDx, 
    const RealVector &v, const RealVector &vDot, PDECoeff &pde) {}
  virtual int getNumODE() const { return 0; }
  struct ODE {
    RealVector c, f;
  };
  virtual void evalODE(double t, const RealVector &v, 
    const RealVector &vdot, 
    const RealMatrix &u, const RealMatrix &DuDx, 
    const RealMatrix &odeR, const RealMatrix &odeDuDt,
    const RealMatrix &odeDuDxDt, RealVector &f) = 0;
  virtual const RealVector &getODEMesh() = 0;

  virtual int getNumEvents() const { return 0; }
  virtual void evalEvents(double t, const RealMatrix &u,
    RealVector &eventsVal, RealVector &eventsIsTerminal,
    RealVector &eventsDirection) = 0;
};

PDESolution pde1d(PDE1dDefn &pde);


