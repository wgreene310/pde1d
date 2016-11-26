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

#include <limits>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

#include <Eigen/LU>

#define DAE_Y_INIT 0

#define SUN_USING_SPARSE 1
#define BAND_SOLVER 0

#include <ida/ida.h>
#if SUN_USING_SPARSE
#include <ida/ida_klu.h>
#endif
#include <ida/ida_band.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

#include "PDE1dImpl.h"
#include "PDE1dOptions.h"
#include "PDE1dException.h"
#include "PDE1dWarningMsg.h"
#include "SunVector.h"
#include "FiniteDiffJacobian.h"
#include "ShapeFunction.h"
#include "ShapeFunctionManager.h"
#include "PDEInitConditions.h"
#include "PDEMeshMapper.h"
#include "PDEModel.h"
#include <util.h>

#define DEBUG_MATS 0

PDE1dImpl::PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options) : 
pde(pde), options(options)
{
#if 0
  sf = std::make_unique<ShapeFunction2>();
#else
  sf = std::unique_ptr<ShapeFunction2>(new ShapeFunction2);
#endif
  sfm = std::unique_ptr<ShapeFunctionManager>(new ShapeFunctionManager);
  ida = 0;
  polyOrder = 1;
  numIntPts = GausLegendreIntRule::getNumPtsForPolyOrder(2 * polyOrder);
  mesh = pde.getMesh();
  checkIncreasing(mesh, 5, "meshPts");
  tspan = pde.getTimeSpan();
  checkIncreasing(tspan, 6, "timePts");
  numNodes = mesh.size();
  numTimes = tspan.size();
  numODE = pde.getNumODE();
  numDepVars = pde.getNumPDE();
  numFEEqns = numNodes*numDepVars;
  totalNumEqns = numFEEqns + numODE;
  dirConsFlagsLeft.resize(numDepVars);
  dirConsFlagsRight.resize(numDepVars);
  y0.resize(totalNumEqns);
  MapMat y0FE(y0.data(), numDepVars, numNodes);
  // get initial conditions
  RealVector ic(numDepVars);
  for (int i = 0; i < numNodes; i++) {
    pde.evalIC(mesh(i), ic);
    y0FE.col(i) = ic;
  }
  pdeModel = std::unique_ptr<PDEModel>(
    new PDEModel(mesh, polyOrder, numDepVars, *sfm.get()));
  if (numODE) {
    v.resize(numODE);
    pde.evalODEIC(v);
    y0.bottomRows(numODE) = v;
    const RealVector &destMesh = pde.getODEMesh();
    meshMapper = std::unique_ptr<PDEMeshMapper>(
    new PDEMeshMapper(mesh, *pdeModel.get(), destMesh));
    size_t numOdePts = destMesh.size();
    vDot.resize(numODE);
    vDot.setZero();
    odeF.resize(numODE);
    odeU.resize(numDepVars, numOdePts);
    odeDuDx.resize(numDepVars, numOdePts);
    odeFlux.resize(numDepVars, numOdePts);
    odeDuDt.resize(numDepVars, numOdePts);
    odeDuDxDt.resize(numDepVars, numOdePts);
  }
  // flag dirichlet constraints
  size_t nm1 = numNodes - 1;
  RealVector ul = y0FE.col(0), ur = y0FE.col(nm1);
  bc.pl.resize(numDepVars);
  bc.pr.resize(numDepVars);
  bc.ql.resize(numDepVars);
  bc.qr.resize(numDepVars);
  pde.evalBC(mesh(0), ul, mesh(nm1), ur, tspan(0), v, vDot, bc);
  for (int i = 0; i < numDepVars; i++) {
    dirConsFlagsLeft[i] = bc.ql[i] == 0;
    dirConsFlagsRight[i] = bc.qr[i] == 0;
  }

  size_t numXPts = 1;
  if (options.isVectorized()) {
    numXPts = numIntPts*(numNodes - 1);
  }
  pdeCoeffs.c.resize(numDepVars, numXPts);
  pdeCoeffs.f.resize(numDepVars, numXPts);
  pdeCoeffs.s.resize(numDepVars, numXPts);

  Cxd.resize(numFEEqns);
  F.resize(numFEEqns);
  S.resize(numFEEqns);

  //printf("Using sparse solver.\n");
  SparseMat P;
  calcJacPattern(P);
  numNonZerosJacMax = P.nonZeros();
  //cout << "P\n" << P << endl;
  finiteDiffJacobian = 
    std::unique_ptr<FiniteDiffJacobian>(new FiniteDiffJacobian(P));

}

PDE1dImpl::~PDE1dImpl()
{
  if(ida) IDAFree(&ida);
}

namespace {

  int resFunc(realtype tres, N_Vector uu, N_Vector up, N_Vector resval,
    void *user_data) {

    PDE1dImpl *pde = (PDE1dImpl*) user_data;

    SunVector u(uu), uDot(up), res(resval);

#define DEBUG 0

#if DEBUG
    printf("resFunc called\n");
    print(uu, "u");
    print(up, "up");
#endif
    pde->calcRHSODE(tres, u, uDot, res);

#if DEBUG
    print(resval, "resval");
    throw PDE1dException("pde1d:debug_exit", "Debugging, quitting early.");
#endif

    return 0;
  }

#if SUN_USING_SPARSE
  int jacFunc(realtype t, realtype c_j,
    N_Vector uu, N_Vector up, N_Vector r,
    SlsMat Jac, void *user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    PDE1dImpl *fem = (PDE1dImpl*) user_data;
    SunVector u(uu), uDot(up), res(r);
    fem->calcJacobianODE(t, c_j, u, uDot, res, Jac);
#if 0
    PrintSparseMat(Jac);
    return 1;
#endif
    return 0;
  }
#endif

  void check_flag(const void *flagvalue, const char *funcname, int opt)
  {

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
      char msg[256];
      sprintf(msg,
        "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      throw PDE1dException("pde1d:sundials_mem_alloc", msg);
    }
    else if (opt == 1) {
      /* Check if flag < 0 */
      int *errflag = (int *)flagvalue;
      if (*errflag < 0) {
        char msg[256];
        sprintf(msg,
          "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
        throw PDE1dException("pde1d:sundials_error", msg);
      }
    }
  }

}

int PDE1dImpl::solveTransient(PDESolution &sol)
{

  sol.time.resize(numTimes);

  const size_t neqImpl = totalNumEqns;
  /* Create vectors uu, up, res, constraints, id. */
  SunVector uu(neqImpl), up(neqImpl), res(neqImpl), id(neqImpl);
  double *udata = &uu[0];
  double *updata = &up[0];
  MapVec u(udata, neqImpl);
  std::copy_n(y0.data(), totalNumEqns, udata);
  up.setConstant(0);
#if 0
  // test for consistent initial conditions
  resFunc(0, uu(), up(), res(), this);
  MapVec resVec(&res[0], neqImpl);
  double initResNorm = resVec.dot(resVec);
  printf("initResNorm=%12.3e\n", sqrt(initResNorm));
#endif
  setAlgVarFlags(up, id);

  /* Call IDACreate and IDAMalloc to initialize solution */
  ida = IDACreate();
  check_flag(ida, "IDACreate", 0);

  int ier = IDASetUserData(ida, this);
  check_flag(&ier, "IDASetUserData", 1);
  ier = IDASetId(ida, id.getNV());
  check_flag(&ier, "IDASetId", 1);
#if 1
  if (numODE) {
    ier = IDASetSuppressAlg(ida, true);
    check_flag(&ier, "IDASetSuppressAlg", 1);
  }
#endif
  double t0 = tspan(0), tf = tspan(numTimes-1);
  PDEInitConditions initCond(ida, *this, uu, up);
  PDEInitConditions::ICPair icPair = initCond.init();
  ier = IDAInit(ida, resFunc, t0, icPair.first->getNV(),
    icPair.second->getNV());
  check_flag(&ier, "IDAInit", 1);
  const double relTol = options.getRelTol(), absTol = options.getAbsTol();
  ier = IDASStolerances(ida, relTol, absTol);
  check_flag(&ier, "IDASStolerances", 1);
  ier = IDASetMaxNumSteps(ida, options.getMaxSteps());
  check_flag(&ier, "IDASetMaxNumSteps", 1);
#if SUN_USING_SPARSE
  //printf("Using sparse solver.\n");
  ier = IDAKLU(ida, (int) totalNumEqns, (int) numNonZerosJacMax);
  check_flag(&ier, "IDAKLU", 1);
  ier = IDASlsSetSparseJacFn(ida, jacFunc);
  check_flag(&ier, "IDASlsSetSparseJacFn", 1);
#elif BAND_SOLVER
  /* Call IDABand to specify the linear solver. */
  int mu = 2*numDepVars-1, ml = mu;
  ier = IDABand(ida, neqImpl, mu, ml);
  check_flag(&ier, "IDABand", 1);
#else
  ier = IDADense(ida, neqImpl);
  check_flag(&ier, "IDADense", 1);
#endif

  // second stage of initial conditions calculation
  initCond.update();
  //initCond.print();

#if 0
  // testing only
  SunVector u0Tmp = initCond.getU0();
  SunVector up0Tmp = initCond.getUp0();
  testICCalc(u0Tmp, up0Tmp, res, id, tf);
#endif
  
  sol.time(0) = tspan(0);
  sol.u.resize(numTimes, numFEEqns);
  sol.u.row(0) = initCond.getU0().topRows(numFEEqns);
  if (numODE) {
    sol.uOde.resize(numTimes, numODE);
    sol.uOde.row(0) = initCond.getU0().bottomRows(numODE);
  }

  for (int i = 1; i < numTimes; i++) {
    double tout=tspan(i), tret;
    ier = IDASolve(ida, tout, &tret, uu.getNV(), up.getNV(), IDA_NORMAL);
    if (ier < 0) {
      printStats();
      char msg[1024];
      sprintf(msg, "Time integration failed at t=%15.6e before reaching final time.\n"
        "Often this is caused by one or more dependent variables becoming unbounded.",
        tret);
      throw PDE1dException("pde1d:integ_failure", msg);
    }
    sol.time(i) = tret;
    sol.u.row(i) = u.topRows(numFEEqns);
    if (numODE)
      sol.uOde.row(i) = u.bottomRows(numODE);
  }
  printStats();

  return 0;
}

void PDE1dImpl::calcGlobalEqns(double time, SunVector &u, SunVector &up,
  RealVector &Cxd, RealVector &F, RealVector &S)
{
  Cxd.setZero(); F.setZero(); S.setZero();
  if (options.isVectorized() && pde.hasVectorPDEEval())
    calcGlobalEqnsVec(time, u, up, Cxd, F, S);
  else
    calcGlobalEqnsScalar(time, u, up, Cxd, F, S);
}

template<class T, class TR>
void PDE1dImpl::calcGlobalEqnsScalar(double t, T &u, T &up,
  TR &Cxd, TR &F, TR &S)
{
  const int nen = numElemNodes; // work around for g++ link error
  const int m = pde.getCoordSystem();
  size_t numElemEqns = numDepVars*nen;
  Cxd.setZero();

  //cout << "u=" << u.transpose() << endl;
  MapMat u2(u.data(), numDepVars, numNodes);
  MapMat up2(up.data(), numDepVars, numNodes);
  //cout << "u2=" << u2.transpose() << endl;
  //cout << "up2=" << up2.transpose() << endl;

  const ShapeFunctionManager::EvaluatedSF &esf = 
    sfm->getShapeFunction(polyOrder);
  const RealMatrix &N = esf.N();
  const RealMatrix &dN = esf.dN();
  const RealVector &intWts = esf.intRuleWts();
#if 0
  cout << "intWts=" << intWts.transpose() << endl;
  cout << "N\n" << N << endl;
#endif

  pdeCoeffs.c.resize(numDepVars, 1);
  pdeCoeffs.f.resize(numDepVars, 1);
  pdeCoeffs.s.resize(numDepVars, 1);

  RealMatrix eCMat = RealMatrix::Zero(numElemEqns, numElemEqns);
  RealVector eF = RealVector::Zero(numElemEqns);
  RealVector eS = RealVector::Zero(numElemEqns);
  RealVector eUp(numElemEqns);
  RealMatrix eCj(nen, nen);
  RealVector eFj(nen), eSj(nen);
  RealVector xiV(1);

  RealVector &x = mesh;
  size_t ne = numNodes - 1;
  size_t eInd = 0;
  for (int e = 0; e < ne; e++) {
    eCMat.setZero();
    eF.setZero();
    eS.setZero();
    double L = x(e + 1) - x(e);
    //fprintf('elem=%d\n', e);
    double jac = L / 2;
    auto dNdx = dN.col(0) / jac;
    for (int i = 0; i < numIntPts; i++) {
      // assume two node elements
      double jacWt = jac*intWts(i);
      double xi = x(e)*N(0, i) + x(e + 1)*N(1, i);
      auto ui = u2.col(e)*N(0, i) + u2.col(e + 1)*N(1, i);
      auto dUiDx = (u2.col(e)*dNdx(0) + u2.col(e + 1)*dNdx(1));
      xiV(0) = xi;
      pde.evalPDE(xiV, t, ui, dUiDx, v, vDot, pdeCoeffs);
      Eigen::Ref<Eigen::VectorXd> cIp = pdeCoeffs.c.col(0);
      Eigen::Ref<Eigen::VectorXd> sIp = pdeCoeffs.s.col(0);
      Eigen::Ref<Eigen::VectorXd> fIp = pdeCoeffs.f.col(0);
      checkCoeffs(pdeCoeffs);
      double xm = 1;
      if (m == 1)
        xm = xi;
      else if(m == 2)
        xm = xi*xi;
      cIp *= xm;
      sIp *= xm;
      fIp *= xm;
 
      for (int j = 0; j < numDepVars; j++) {
        eCj  = N.col(i)*cIp(j)*N.col(i).transpose()*jacWt;
        //eCMat += eCj;
        eFj = dNdx*fIp(j)*jacWt;
        //eF += eFj;
        eSj = N.col(i)*sIp(j)*jacWt;
        //eS += eSj;
        for (int k = 0; k < numElemNodes; k++) {
          size_t kk = j + k*numDepVars;
          eF(kk) += eFj(k);
          eS(kk) += eSj(k);
          for (int l = 0; l < numElemNodes; l++) {
            size_t ll = j + l*numDepVars;
            eCMat(kk, ll) += eCj(k, l);
          }
        }
      }
    }
#if DEBUG_MATS
    cout << "eF: " << eF.transpose() << endl;
    cout << "eCMat\n" << eCMat << endl;
#endif
    eUp.head(numDepVars) = up2.col(e);
    eUp.tail(numDepVars) = up2.col(e+1);
#if DEBUG_MATS
    //cout << "eUp=" << eUp.transpose() << endl;
#endif
    Cxd.segment(eInd, numElemEqns) += eCMat*eUp;
    F.segment(eInd, numElemEqns) += eF;
    S.segment(eInd, numElemEqns) += eS;
    eInd += numDepVars;
  }
  //cout << "F\n" << F << endl;
}

  template<class T, class TR>
  void PDE1dImpl::calcGlobalEqnsVec(double t, T &u, T &up,
    TR &Cxd, TR &F, TR &S)
  {
    const int nen = numElemNodes; // work around for g++ link error
    const int m = pde.getCoordSystem();
    size_t numElemEqns = numDepVars*nen;
    Cxd.setZero();

    //cout << "u=" << u.transpose() << endl;
    MapMat u2(u.data(), numDepVars, numNodes);
    MapMat up2(up.data(), numDepVars, numNodes);
    //cout << "u2=" << u2.transpose() << endl;
    //cout << "up2=" << up2.transpose() << endl;

    const ShapeFunctionManager::EvaluatedSF &esf =
      sfm->getShapeFunction(1);
    const RealMatrix &N = esf.N();
    const RealMatrix &dN = esf.dN();
    const RealVector &intWts = esf.intRuleWts();
#if 0
    cout << "intWts=" << intWts.transpose() << endl;
    cout << "N\n" << N << endl;
#endif

    // get the coefficients at all integ pts in a single call
    size_t ne = numNodes - 1;
    size_t numXPts = numIntPts*ne;
    xPts.resize(numXPts);
    uPts.resize(numDepVars, numXPts);
    duPts.resize(numDepVars, numXPts);

    RealVector &x = mesh;
    int ip = 0;
    for (int e = 0; e < ne; e++) {
      double L = x(e + 1) - x(e);
      double jac = L / 2;
      for (int i = 0; i < numIntPts; i++) {
        // assume two node elements
        double jacWt = jac*intWts(i);
        xPts(ip) = x(e)*N(0, i) + x(e + 1)*N(1, i);
        uPts.col(ip) = u2.col(e)*N(0, i) + u2.col(e + 1)*N(1, i);
        auto dNdx = dN.col(0) / jac;
        duPts.col(ip) = (u2.col(e)*dNdx(0) + u2.col(e + 1)*dNdx(1));
        ip++;
      }
    }
    pdeCoeffs.c.resize(numDepVars, numXPts);
    pdeCoeffs.f.resize(numDepVars, numXPts);
    pdeCoeffs.s.resize(numDepVars, numXPts);
    pde.evalPDE(xPts, t, uPts, duPts, v, vDot, pdeCoeffs);

    RealMatrix eCMat = RealMatrix::Zero(numElemEqns, numElemEqns);
    RealVector eF = RealVector::Zero(numElemEqns);
    RealVector eS = RealVector::Zero(numElemEqns);
    RealVector eUp(numElemEqns);
    RealMatrix eCj(nen, nen);
    RealVector eFj(nen), eSj(nen);

    ip = 0; 
    size_t eInd = 0;
    for (int e = 0; e < ne; e++) {
      eCMat.setZero();
      eF.setZero();
      eS.setZero();
      double L = x(e + 1) - x(e);
      //fprintf('elem=%d\n', e);
      double jac = L / 2;
      auto dNdx = dN.col(0) / jac;
      for (int i = 0; i < numIntPts; i++) {
        // assume two node elements
        double jacWt = jac*intWts(i);
        double xi = xPts(ip);
        Eigen::Ref<Eigen::VectorXd> cIp = pdeCoeffs.c.col(ip);
        Eigen::Ref<Eigen::VectorXd> sIp = pdeCoeffs.s.col(ip);
        Eigen::Ref<Eigen::VectorXd> fIp = pdeCoeffs.f.col(ip);
        double xm = 1;
        if (m == 1)
          xm = xi;
        else if (m == 2)
          xm = xi*xi;
        cIp *= xm;
        sIp *= xm;
        fIp *= xm;

        for (int j = 0; j < numDepVars; j++) {
          eCj = N.col(i)*cIp(j)*N.col(i).transpose()*jacWt;
          //eCMat += eCj;
          eFj = dNdx*fIp(j)*jacWt;
          //eF += eFj;
          eSj = N.col(i)*sIp(j)*jacWt;
          //eS += eSj;
          for (int k = 0; k < numElemNodes; k++) {
            size_t kk = j + k*numDepVars;
            eF(kk) += eFj(k);
            eS(kk) += eSj(k);
            for (int l = 0; l < numElemNodes; l++) {
              size_t ll = j + l*numDepVars;
              eCMat(kk, ll) += eCj(k, l);
            }
          }
        }
        ip++;
      }
#if DEBUG_MATS
      cout << "eF: " << eF.transpose() << endl;
      cout << "eCMat\n" << eCMat << endl;
#endif
      eUp.head(numDepVars) = up2.col(e);
      eUp.tail(numDepVars) = up2.col(e + 1);
#if DEBUG_MATS
      //cout << "eUp=" << eUp.transpose() << endl;
#endif
      Cxd.segment(eInd, numElemEqns) += eCMat*eUp;
      F.segment(eInd, numElemEqns) += eF;
      S.segment(eInd, numElemEqns) += eS;
      eInd += numDepVars;
    }
    //cout << "F\n" << F << endl;
  }

  void PDE1dImpl::calcRHSODE(double time, SunVector &u, SunVector &up, 
    SunVector &R)
{
    // copy the ode dofs to their own vectors
    if (numODE) {
      v = u.bottomRows(numODE);
      vDot = up.bottomRows(numODE);
    }

    calcGlobalEqns(time, u, up, Cxd, F, S);
    R.topRows(numFEEqns) = F - S;

  // apply constraints
  MapMat u2(u.data(), numDepVars, numNodes);
  size_t rightDofOff = numFEEqns - numDepVars;
  size_t nm1 = numNodes - 1;
  RealVector ul = u2.col(0), ur = u2.col(nm1);
  pde.evalBC(mesh(0), ul, mesh(nm1), ur, time, v, vDot, bc);
  for (int i = 0; i < numDepVars; i++) {
    if (bc.ql(i) != 0) {
      R(i) -= bc.pl(i) / bc.ql(i);
    }
    else {
      // apply dirichlet constraint
      R(i) = bc.pl(i);
      Cxd(i) = 0;
    }
    if (bc.qr(i) != 0) {
      R(i + rightDofOff) += bc.pr(i) / bc.qr(i);
    }
    else {
      // apply dirichlet constraint
      R(i + rightDofOff) = bc.pr(i);
      Cxd(i + rightDofOff) = 0;
    }
  }
  R.topRows(numFEEqns) += Cxd;
#if 0
  cout << "Cxd\n" << Cxd << endl;
  cout << "F\n" << F << endl;
#endif
  // add odes, if any
  if (numODE) {
    MapMat f2(F.data(), numDepVars, numNodes);
    meshMapper->mapFunction(u2, odeU);
    meshMapper->mapFunctionDer(u2, odeDuDx);
    meshMapper->mapFunctionDer(f2, odeFlux);
    MapMat up2(up.data(), numDepVars, numNodes);
    meshMapper->mapFunction(up2, odeDuDt);
    meshMapper->mapFunctionDer(up2, odeDuDxDt);
    pde.evalODE(time, v, vDot, odeU, odeDuDx, odeFlux, 
      odeDuDt, odeDuDxDt, odeF);
    R.bottomRows(numODE) = odeF;
  }
}

void PDE1dImpl::testMats()
{
  SunVector u(numFEEqns), up(numFEEqns);
  for (int i = 0; i < numFEEqns; i++) {
    u(i) = .1*(i+1);
    up(i) = 0;
  }
  RealVector Cxd = RealVector::Zero(numFEEqns);
  RealVector F = RealVector::Zero(numFEEqns);
  RealVector S = RealVector::Zero(numFEEqns);
  calcGlobalEqns(0, u, up, Cxd, F, S);
  cout << "Cxd\n" << Cxd.transpose() << endl;
  cout << "F\n" << F.transpose() << endl;
  cout << "S\n" << S.transpose() << endl;
#if 0
  RealVector R = RealVector::Zero(numFEMEqns);
  calcRHSODE(0, u, up, R);
#else
  RealVector R = S - F;
#endif
  cout << "R\n" << R.transpose() << endl;
  throw PDE1dException("pde1d:testResidual", "residual test complete.");
}


void PDE1dImpl::setAlgVarFlags(SunVector &y0p, SunVector &id)
{

  MapMat y0FE(y0.data(), numDepVars, numNodes);
  
  double t0 = tspan[0];

  id.setConstant(1);
  double *iddata = id.data();
  MapMat idMat(iddata, numDepVars, numNodes);

  // flag equations where c==0
  const int nen = numElemNodes;
  RealMatrix dN(nen, nen);
  double xi[] = { -1, 1 };
  for (int i = 0; i < numElemNodes; i++) {
    double pt = xi[i];
    sf->dNdr(pt, dN.col(i).data());
  }
 
  duPts.resize(numDepVars, numNodes);
  size_t nnm1 = numNodes - 1;
  for(int i = 0; i < numNodes; i++) {
    double L;
    if (i < nnm1)
      L = mesh(i + 1) - mesh(i);
    else
      L = mesh(i) - mesh(i - 1);
    double jac = L / 2;
    auto dNdx = dN.col(0) / jac;
    if (i < nnm1)
      duPts.col(i) = (y0FE.col(i)*dNdx(0) + y0FE.col(i + 1)*dNdx(1));
    else
      duPts.col(i) = (y0FE.col(i-1)*dNdx(0) + y0FE.col(i)*dNdx(1));
  }

  size_t numXPts = 1;
  if (options.isVectorized())
    numXPts = numNodes;
  pdeCoeffs.c.resize(numDepVars, numXPts);
  pdeCoeffs.f.resize(numDepVars, numXPts);
  pdeCoeffs.s.resize(numDepVars, numXPts);

  if (options.isVectorized() && pde.hasVectorPDEEval()) {
    pde.evalPDE(mesh, 0, y0FE, duPts, 
      v, vDot, pdeCoeffs);
    for (int i = 0; i < numNodes; i++) {
      for (int j = 0; j < numDepVars; j++)
        if (pdeCoeffs.c(j,i)==0)
          idMat(j, i) = 0;
    }
  }
  else {
    RealVector xiV(1);
    for (int i = 0; i < numNodes; i++) {
      double xi = mesh(i);
      const auto &ui = y0FE.col(i);
      const auto &dUiDx = duPts.col(i);
      xiV(0) = xi;
      Eigen::Ref<Eigen::VectorXd> cIp = pdeCoeffs.c.col(0);
      pde.evalPDE(xiV, t0, ui, dUiDx, v, vDot, pdeCoeffs);
      for (int j = 0; j < numDepVars; j++)
        if (cIp[j] == 0)
          idMat(j, i) = 0;
    }
  }

  // account for dirichlet constraints at ends
  size_t rtBcOff = numFEEqns - numDepVars;
  for (int i = 0; i < numDepVars; i++) {
    if (dirConsFlagsLeft[i])
      iddata[i] = 0;
    if (dirConsFlagsRight[i])
      iddata[i + rtBcOff] = 0;
  }

  // check ODEs
  if (numODE) {
    SunVector y0Tmp(totalNumEqns);
    y0Tmp = y0;
    calcGlobalEqns(t0, y0Tmp, y0p, Cxd, F, S);
    MapMat  f2(F.data(), numDepVars, numNodes);
    MapMat  yp0FE(y0p.data(), numDepVars, numNodes);
    v = y0.bottomRows(numODE);
    vDot = y0p.bottomRows(numODE);
    RealMatrix odeJac = calcODEJacobian(t0, y0FE, yp0FE, f2, v, vDot);
    //cout << "odeJac:\n" << odeJac << endl;
    for (int i = 0; i < numODE; i++) {
      auto rci = odeJac.row(i) + odeJac.col(i).transpose();
      if ((rci.array() == 0).all())
        id(numFEEqns + i) = 0;
    }
    //print(id.getNV(), "id");
  }
}

RealMatrix PDE1dImpl::calcODEJacobian(double time, const RealMatrix &yFE, 
  const RealMatrix &ypFE, const RealMatrix &f2, RealVector &v, RealVector &vdot)
{
  RealMatrix jac(numODE, numODE);
  meshMapper->mapFunction(yFE, odeU);
  meshMapper->mapFunctionDer(yFE, odeDuDx);
  meshMapper->mapFunctionDer(f2, odeFlux);

  meshMapper->mapFunction(ypFE, odeDuDt);
  meshMapper->mapFunctionDer(ypFE, odeDuDxDt);
  pde.evalODE(time, v, vDot, odeU, odeDuDx, odeFlux, 
    odeDuDt, odeDuDxDt, odeF);
  RealVector vdotDelta = vDot, fDelta(numODE);
  double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
  for (int i = 0; i < numODE; i++) {
    double h = sqrtEps*std::max(vDot[i], 1.0);
    double tmpVdotI = vDot[i];
    vdotDelta[i] += h;
    pde.evalODE(time, v, vdotDelta, odeU, odeDuDx, odeFlux, 
      odeDuDt, odeDuDxDt, fDelta);
    jac.col(i) = (fDelta - odeF) / h;
    vdotDelta[i] = tmpVdotI;
  }
  return jac;
}

void PDE1dImpl::checkIncreasing(const RealVector &v, 
  int argNum, const char *argName)
{
  for (int i = 1; i < v.size(); i++) {
    if (v[i] <= v[i - 1]) {
      char msg[1024];
      sprintf(msg, "The values in argument number %d (\"%s\") "
        "must be strictly increasing.\n"
        "Entry %d is %g but entry %d is %g.",
        argNum, argName, i, v[i - 1], i + 1, v[i]);
      throw PDE1dException("pde1d:decreasing_indep_var", msg);
    }
  }
}

void PDE1dImpl::checkCoeffs(const PDE1dDefn::PDECoeff &coeffs)
{
  for (int j = 0; j < coeffs.c.cols(); j++) {
    int zeroCol = -1;
    bool allZero = true;
    for (int i = 0; i < coeffs.c.rows(); i++) {
      if (coeffs.c.col(j)[i] != 0) {
        allZero = false;
        break;
      }
    }
    if (allZero) {
      throw PDE1dException("pde1d:no_parabolic_eqn", "At least one of the "
        "entries in the c-coefficient vector must be non-zero.");
    }
  }
}

void PDE1dImpl::printStats()
{
  if (!options.printStats()) return;
  long nsteps, nrevals, nlinsetups, netfails;
  int klast, kcur;
  double hinused, hlast, hcur, tcur;
  int flag = IDAGetIntegratorStats(ida, &nsteps, &nrevals, &nlinsetups,
    &netfails, &klast, &kcur, &hinused,
    &hlast, &hcur, &tcur);
  long nniters, nncfails;
  flag = IDAGetNonlinSolvStats(ida, &nniters, &nncfails);
  printf("Number of internal time steps = %ld\n", nsteps);
  printf("Number of residual function calls = %ld\n", nrevals);
  printf("Number of Jacobian calculations = %ld\n", nlinsetups);
  printf("Number of solution accuracy test failures = %ld\n", netfails);
  printf("Last internal time step size = %12.3e\n", hlast);
  printf("Number of nonlinear iterations = %ld\n", nniters);
  printf("Number of nonlinear convergence failures = %ld\n", nncfails);
}

void PDE1dImpl::calcJacPattern(Eigen::SparseMatrix<double> &J)
{
  J.resize(totalNumEqns, totalNumEqns);
  size_t n2 = numDepVars*numDepVars;
  size_t nel = numNodes - 1;
  size_t nnz = 3 * n2*(nel + 1); // approximate nnz
  J.reserve(nnz);
  size_t eOff = 0;
  for (int n = 0; n < nel + 1; n++) {
    for (int i = 0; i < numDepVars; i++) {
      for (int j = 0; j < numDepVars; j++) {
        J.insert(i + eOff, j + eOff) = 1;
        if (n < nel) {
          J.insert(i + eOff + numDepVars, j + eOff) = 1;
          J.insert(i + eOff, j + eOff + numDepVars) = 1;
        }
      }
    }
    eOff += numDepVars;
  }
  // odes may couple with any other vars
  eOff = numNodes*numDepVars;
  for (int i = 0; i < numODE; i++) {
    for (int j = 0; j < numODE; j++)
      J.insert(i + eOff, j + eOff) = 1;
    for (int j = 0; j < numFEEqns; j++) {
      J.insert(i + eOff, j) = 1;
      J.insert(j, i + eOff) = 1;
    }
  }
  J.makeCompressed();
  //cout << "pattern\n" << J.toDense() << endl;
}

void PDE1dImpl::calcJacobianODE(double time, double beta, SunVector &u, 
  SunVector &up, SunVector &res, SlsMat Jac)
{
  finiteDiffJacobian->calcJacobian(time, 1, beta, u.getNV(), up.getNV(), 
    res.getNV(), resFunc, this, Jac);
}

void PDE1dImpl::calcJacobian(double time, double alpha, double beta, SunVector &u,
  SunVector &up, SunVector &R, SparseMat &Jac)
{
  const bool useCD = !true;
  finiteDiffJacobian->calcJacobian(time, alpha, beta, u.getNV(), up.getNV(), 
    R.getNV(), resFunc, this, Jac, useCD);
}

void PDE1dImpl::testICCalc(SunVector &uu, SunVector &up, SunVector &res,
  SunVector &id, double tf)
{
  resFunc(0, uu.getNV(), up.getNV(), res.getNV(), this);
  print(uu.getNV(), "u");
  print(up.getNV(), "up");
  print(res.getNV(), "res");
  // calc jac matrices
  SparseMat dfDy(totalNumEqns, totalNumEqns), 
    dfDyp(totalNumEqns, totalNumEqns);
  finiteDiffJacobian->calcJacobian(0, 1, 0, uu.getNV(), up.getNV(), 
    res.getNV(), resFunc, this, dfDy);
  cout << "dfDy\n" << dfDy.toDense() << endl;
  Eigen::FullPivLU<RealMatrix> lu(dfDy.toDense());
  printf("dfdy: n=%d, rank=%d\n", dfDy.rows(), lu.rank());

#if 0
  fDiffJac->calcJacobian(0, 0, 1, uu.getNV(), up.getNV(), res.getNV(), 
    resFunc, this, dfDyp);
#else
  calcJacobian(0, 0, 1, uu, up, res, dfDyp);
#endif
  //cout << "dfDyp\n" << dfDyp.toDense() << endl;
  printMat(dfDyp.toDense(), "dfDyp");
#if 0

  int ier = IDACalcIC(ida, IDA_YA_YDP_INIT, tf);
  print(id.getNV(), "id vector");
  SunVector yy0_mod(totalNumEqns), yp0_mod(totalNumEqns);
  ier = IDAGetConsistentIC(ida, yy0_mod.getNV(), yp0_mod.getNV());
  print(yy0_mod.getNV(), "ICy");
  print(yp0_mod.getNV(), "ICyp");
#endif

  throw PDE1dException("pde1d:test", "IC test completed");
}

double PDE1dImpl::calcResidualNorm(double t, SunVector &uu, SunVector &up, SunVector &res)
{
  resFunc(t, uu.getNV(), up.getNV(), res.getNV(), this);
  MapVec r(&res[0], res.size());
  double rRms = sqrt(r.dot(r));
  return rRms;
}