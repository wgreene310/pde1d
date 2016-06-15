// Copyright (C) 2016 William H. Greene
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

#ifndef PDE1dImpl_h
#define PDE1dImpl_h

#include <vector>

#include "PDE1dDefn.h"
#include "GausLegendreIntRule.h"

#include <nvector/nvector_serial.h>

class PDE1dOptions;

class PDE1dImpl {
public:
  PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options);
  ~PDE1dImpl();
  int solveTransient(PDESolution &sol);
  template<class T>
  void calcRHSODE(double time, T &u, T &up, T &R);
  void testMats();
private:
  template<class T, class TR>
  void calcGlobalEqns(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  void setAlgVarFlags(N_Vector id);
  void checkIncreasing(const RealVector &v, int argNum, const char *argName);
  PDE1dDefn &pde;
  PDE1dOptions &options;
  GausLegendreIntRule *intRule;
  RealVector mesh, tspan;
  int numNodes, numTimes;
  int numDepVars, numFEMEqns;
  static const int numElemNodes = 2;
  std::vector<bool> dirConsFlagsLeft, dirConsFlagsRight;
  RealMatrix y0;
  PDE1dDefn::BC bc;
  PDE1dDefn::PDE coeffs;
};

#define FUNC_NAME "pde1d"

#endif

