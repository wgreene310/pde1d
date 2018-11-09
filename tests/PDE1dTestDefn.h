#pragma once
#include "PDE1dDefn.h"

class PDE1dTestDefn : public PDE1dDefn
{
public:
  PDE1dTestDefn(double L, int nel, double tFinal,
    int nt);
  virtual const RealVector &getMesh() const {
    return mesh;
  }
  virtual const RealVector &getTimeSpan() const {
    return tspan;
  }
  virtual void evalODEIC(RealVector &ic) { }
  virtual const RealVector &getODEMesh() {
    return odeMesh;
  }
private:
  RealVector mesh, tspan, odeMesh;
};