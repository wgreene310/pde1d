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
class ShapeFunctionManager;
class SunVector;
struct _SlsMat;
class PDEMeshMapper;
class PDEModel;

class PDE1dImpl {
public:
  PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options);
  ~PDE1dImpl();
  int solveTransient(PDESolution &sol);
  void calcRHSODE(double time, SunVector &u, SunVector &up, SunVector &R);
  void calcJacobianODE(double time, double alpha, SunVector &u, 
    SunVector &up, SunVector &R, _SlsMat *Jac);
  void calcJacobian(double time, double alpha, double beta, SunVector &u,
    SunVector &up, SunVector &R, SparseMat &Jac);
  const PDE1dOptions &getOptions() const {
    return options;
  }
  const RealVector &getTimePoints() const {
    return tspan;
  }
  void testMats();
private:
  void calcGlobalEqns(double t, SunVector &u, SunVector &up, 
    RealVector &Cxd, RealVector &F, RealVector &S);
  template<class T, class TR>
  void calcGlobalEqnsScalar(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  template<class T, class TR>
  void calcGlobalEqnsVec(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  void setAlgVarFlags(SunVector &y0p, SunVector &id);
  RealMatrix calcODEJacobian(double time, const RealMatrix &yFE, 
    const RealMatrix &ypFE, const RealMatrix &r2, RealVector &v, RealVector &vdot);
  void checkIncreasing(const RealVector &v, int argNum, const char *argName);
  void checkCoeffs(const PDE1dDefn::PDECoeff &coeffs);
  void printStats();
  void calcJacPattern(Eigen::SparseMatrix<double> &jac);
  void testICCalc(SunVector &uu, SunVector &up, SunVector &res,
    SunVector &id, double tf);
  double calcResidualNorm(double t, SunVector &uu, SunVector &up, SunVector &res);
  PDE1dDefn &pde;
  PDE1dOptions &options;
  RealVector mesh, tspan;
  size_t numNodes, numTimes;
  size_t numDepVars, numODE, numFEEqns, totalNumEqns;
  size_t numNonZerosJacMax;
  static const int numElemNodes = 2;
  int polyOrder, numIntPts;
  std::vector<bool> dirConsFlagsLeft, dirConsFlagsRight;
  RealVector y0;
  PDE1dDefn::BC bc;
  PDE1dDefn::PDECoeff pdeCoeffs;
  RealVector Cxd, F, S;
  // temporary arrays for vectorized mode
  RealVector xPts;
  RealMatrix uPts, duPts;
  std::unique_ptr<FiniteDiffJacobian > finiteDiffJacobian;
  void *ida;
  std::unique_ptr<ShapeFunction> sf;
  std::unique_ptr<ShapeFunctionManager> sfm;
  std::unique_ptr<PDEMeshMapper> meshMapper;
  std::unique_ptr<PDEModel> pdeModel;
  RealVector v, vDot, odeF;
  RealMatrix odeU, odeDuDx, odeFlux, odeDuDt, odeDuDxDt;
};

#define FUNC_NAME "pde1d"

#endif

