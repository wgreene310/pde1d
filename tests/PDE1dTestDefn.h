#pragma once
#include "PDE1dDefn.h"

class PDE1dTestDefn : public PDE1dDefn
{
public:
  PDE1dTestDefn(double L, int nel, double tFinal,
    int nt, int numPde=1, int numOde=0, int numOdeXPts=0);
  virtual const RealVector &getMesh() const {
    return mesh;
  }
  virtual const RealVector &getTimeSpan() const {
    return tspan;
  }
  virtual const RealVector &getODEMesh() {
    return odeMesh;
  }
  virtual void evalIC(double x, RealVector &ic) {
    ic.setZero();
  }
  virtual void evalODEIC(RealVector &ic) {
    ic.setZero();
  }
  virtual int getNumPDE() const { return numPde; }
  virtual int getNumODE() const { return numOde; }
protected:
  RealVector mesh, tspan, odeMesh;
  int numPde, numOde;
};