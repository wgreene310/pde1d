// Copyright (C) 2016-2017 William H. Greene
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
#if SUNDIALS_3
#include <sundials/sundials_matrix.h>
#else
#include <sundials/sundials_sparse.h>
#endif

class PDE1dOptions;
class FiniteDiffJacobian;
class ShapeFunction;
class ShapeFunctionManager;
class SunVector;
class PDEMeshMapper;
class PDEModel;
class PDEEvents;

class PDE1dImpl {
public:
  PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options);
  ~PDE1dImpl();
  int solveTransient(PDESolution &sol);
  void calcRHSODE(double time, SunVector &u, SunVector &up, SunVector &R);
#if SUNDIALS_3
  void calcJacobianODE(double time, double alpha, SunVector &u, 
    SunVector &up, SunVector &R, SUNMatrix Jac);
#else
void calcJacobianODE(double time, double alpha, SunVector &u, 
    SunVector &up, SunVector &R, SlsMat Jac);
#endif
  void calcJacobian(double time, double alpha, double beta, SunVector &u,
    SunVector &up, SunVector &R, SparseMat &Jac);
  void calcEvents(double time, const SunVector &u, double *gOut);
  const PDE1dOptions &getOptions() const {
    return options;
  }
  const RealVector &getTimePoints() const {
    return tspan;
  }
  const PDEModel &getModel() const {
    return *pdeModel;
  }
private:
  void calcGlobalEqns(double t, SunVector &u, SunVector &up, 
    RealVector &Cxd, RealVector &F, RealVector &S);
  template<class T, class TR>
  void calcGlobalEqnsNonVectorized(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  template<class T, class TR>
  void calcGlobalEqnsVectorized(double t, T &u, T &up, TR &Cxd, TR &F, TR &S);
  void setAlgVarFlags(SunVector &y0, SunVector &y0p, SunVector &id);
  RealMatrix calcDOdeDvDot(double time, const RealMatrix &yFE, 
    const RealMatrix &ypFE, const RealMatrix &r2, RealVector &v, RealVector &vdot);
  RealMatrix calcDOdeDv(double time, const RealMatrix &yFE,
    const RealMatrix &ypFE, const RealMatrix &r2, RealVector &v, RealVector &vdot);
  void calcDOdeDu(double time, const RealMatrix &yFE,
    const RealMatrix &ypFE, const RealMatrix &r2, RealVector &v, 
    RealVector &vdot, RealMatrix &jac, RealMatrix &jacDot);
  void checkIncreasing(const RealVector &v, int argNum, const char *argName);
  void checkCoeffs(const PDE1dDefn::PDECoeff &coeffs);
  void printStats();
  void calcJacPattern(Eigen::SparseMatrix<double> &jac);
  void testICCalc(SunVector &uu, SunVector &up, SunVector &res,
    SunVector &id, double tf);
  double calcResidualNorm(double t, SunVector &uu, SunVector &up, SunVector &res);
  void getFEInitConditions(RealVector &y0);
  template<class T>
  void interpolateGlobalVecToViewMesh(const T &gVec,
    RealMatrix &viewVec);
  void jacobianDiagnostics(double t0, SunVector &u,
     SunVector &up, SunVector &R);
  void testMats(const RealVector &y0);
  void testODEJacobian(RealVector &y);
  PDE1dDefn &pde;
  PDE1dOptions &options;
  RealVector mesh, tspan;
  size_t numTimes;
  size_t numDepVars, numODE, numFEEqns, totalNumEqns;
  size_t numNonZerosJacMax;
  int polyOrder, numIntPts;
  std::vector<bool> dirConsFlagsLeft, dirConsFlagsRight;
  PDE1dDefn::BC bc;
  PDE1dDefn::PDECoeff pdeCoeffs;
  RealVector Cxd, F, S;
  // temporary arrays for vectorized mode
  RealVector xPts;
  RealMatrix uPts, duPts;
  std::unique_ptr<FiniteDiffJacobian> finiteDiffJacobian;
  void *ida;
  std::unique_ptr<ShapeFunction> sf;
  std::unique_ptr<ShapeFunctionManager> sfm;
  std::unique_ptr<PDEMeshMapper> meshMapper;
  std::unique_ptr<PDEModel> pdeModel;
  std::unique_ptr<PDEEvents> pdeEvents;
  RealVector v, vDot, odeF;
  RealMatrix odeU, odeDuDx, odeFlux, odeDuDt, odeDuDxDt;
  IntVector isOdeAConstraint;
  size_t numViewElemsPerElem;
};

#define FUNC_NAME "pde1d"

#endif

