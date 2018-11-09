
#include "PDE1dTestDefn.h"

PDE1dTestDefn::PDE1dTestDefn(double L, int nel, double tFinal, int nt)
{
  mesh.resize(nel + 1);
  tspan.resize(nt);
  double dx = L / nel, x = 0;
  for (int i = 0; i < nel + 1; i++) {
    mesh[i] = x;
    x += dx;
  }
  double dt = tFinal / (nt - 1), t = 0;
  for (int i = 0; i < nt; i++) {
    tspan[i] = t;
    t += dt;
  }
}