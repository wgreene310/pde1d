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
#include <memory>

#include <Eigen/SparseCore>
typedef Eigen::SparseMatrix<double> SparseMat;

#include "PDE1dDefn.h"
#include "GausLegendreIntRule.h"

#include <nvector/nvector_serial.h>

class PDE1dOptions;
class FiniteDiffJacobian;
class ShapeFunction;

class PDE1dImpl {
public:
  PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options);
  ~PDE1dImpl();
  int solveTransient(PDESolution &sol);
  template<class T>
  void calcRHSODE(double time, T &u, T &up, T &R);
  template<class T, class T2>
  void calcJacobianODE(double time, double alpha, T &u, T &up, T &R,
    T2 Jac);
  void testMats();
private:
  template<class T, class TR>
  void calcGlobalEqns(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  template<class T, class TR>
  void calcGlobalEqnsVec(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  void setAlgVarFlags(N_Vector id);
  void checkIncreasing(const RealVector &v, int argNum, const char *argName);
  void checkCoeffs(const PDE1dDefn::PDE &coeffs);
  void printStats();
  void calcJacPattern(Eigen::SparseMatrix<double> &jac);
  PDE1dDefn &pde;
  PDE1dOptions &options;
  GausLegendreIntRule *intRule;
  RealVector mesh, tspan;
  int numNodes, numTimes;
  int numDepVars, numODE, numFEMEqns;
  static const int numElemNodes = 2;
  std::vector<bool> dirConsFlagsLeft, dirConsFlagsRight;
  RealMatrix y0;
  PDE1dDefn::BC bc;
  PDE1dDefn::PDE coeffs;
  PDE1dDefn::PDEVec coeffsAllPts; // FIXME should have only one of coeffs or this one
  // temporary arrays for vectorized mode
  RealVector xPts;
  RealMatrix uPts, duPts;
  FiniteDiffJacobian *fDiffJac;
  void *ida;
  std::unique_ptr<ShapeFunction> sf;
};

#define FUNC_NAME "pde1d"

#endif

