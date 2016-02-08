#pragma once

#include <Eigen/Core>

typedef Eigen::VectorXd RealVector;
typedef Eigen::MatrixXd RealMatrix;

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
};

struct PDESolution {
  PDESolution(int nx=1, int nt=1) : u(nt, nx), time(nt) {}
  RealMatrix u;
  RealVector time;
};

PDESolution pde1d(PDE1dDefn &pde);

