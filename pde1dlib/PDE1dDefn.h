#pragma once

#include <MatrixTypes.h>

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
  virtual RealVector getMesh() = 0;
  virtual RealVector getTimeSpan() = 0;
  virtual bool hasVectorPDEEval() const { return false;  }
  struct PDEVec {
    RealMatrix c, f, s;
  };
  virtual void evalPDE(const RealVector &x, double t,
    const RealMatrix &u, const RealMatrix &DuDx, 
    const RealVector &v, const RealVector &vDot, PDEVec &pde) {}
  virtual int getNumODE() const { return 0; }
  struct ODE {
    RealVector c, f;
  };
  virtual void evalODE(double t, const RealVector &v, 
    const RealVector &vdot, 
    const RealMatrix &u, const RealMatrix &DuDx, RealVector &f) = 0;
  virtual const RealVector &getODEMesh() = 0;
};

struct PDESolution {
  PDESolution(int nx=1, int nt=1) {}
  RealMatrix u;
  RealMatrix uOde;
  RealVector time;
};

PDESolution pde1d(PDE1dDefn &pde);

