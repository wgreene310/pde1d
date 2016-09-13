#pragma once

#include <MatrixTypes.h>

class PDE1dDefn
{
public:
  PDE1dDefn();
  virtual int getNumEquations() = 0;
  virtual int getCoordSystem() { return 0; }
  virtual void evalIC(double x, RealVector &ic) = 0;
  struct BC {
    RealVector pl, ql, pr, qr;
  };
  virtual void evalBC(double xl, const RealVector &ul,
    double xr, const RealVector &ur, double t, BC &bc) = 0;
  struct PDE {
    RealVector c, f, s;
  };
  virtual void evalPDE(double x, double t, 
    const RealVector &u, const RealVector &DuDx, PDE &pde) = 0;
  virtual RealVector getMesh() = 0;
  virtual RealVector getTimeSpan() = 0;
  virtual bool hasVectorPDEEval() const { return false;  }
  struct PDEVec {
    RealMatrix c, f, s;
  };
  virtual void evalPDE(RealVector x, double t,
    const RealMatrix &u, const RealMatrix &DuDx, PDEVec &pde) {}
  virtual bool hasODE() const { return false; }
  virtual int getNumODE() const { return 0; }
  struct ODE {
    RealVector c, f;
  };
};

struct PDESolution {
  PDESolution(int nx=1, int nt=1) : u(nt, nx), time(nt) {}
  RealMatrix u;
  RealVector time;
};

PDESolution pde1d(PDE1dDefn &pde);

