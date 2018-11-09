#include <iostream>

#pragma once

#include "PDE1dTestDefn.h"

class ExampleHeatCond : public PDE1dTestDefn
{
public:
  ExampleHeatCond(double L, int nel, double tFinal,
    int nt) :  L(L), nel(nel), tFinal(tFinal), nt(nt),
    PDE1dTestDefn(L, nel, tFinal, nt) { }
  virtual int getNumPDE() const { return 1;  }
  virtual void evalIC(double x, RealVector &ic) {
    ic(0) = 0;
  }
  virtual void evalBC(double xl, const RealVector &ul,
    double xr, const RealVector &ur, double t, 
    const RealVector &v, const RealVector &vDot, BC &bc) {
    bc.pl << ul(0) - 10;
    bc.ql << 0;
#if 1
    // insulated at rt end
    bc.pr << 0;
    bc.qr << 1;
#else
    // dirichlet at rt end
    bc.pr << ur(0) - 20;
    bc.qr << 0;
#endif
  }
  virtual void evalPDE(double x, double t,
    const RealVector &u, const RealVector &DuDx, 
    const RealVector &v, const RealVector &vDot, PDECoeff &pde) {
    pde.c(0) = 1;
    pde.f = 10 * DuDx;
    pde.s(0) = 0;
  }
private:
  double L, tFinal;
  int nel, nt;
  //RealVector mesh, tspan, odeMesh;
};

