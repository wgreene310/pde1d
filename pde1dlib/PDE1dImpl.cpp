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

#include <limits>
#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
#include <algorithm>
#include <memory>

#include <Eigen/LU>

#define DEBUG_MATS 0
#define DAE_Y_INIT 0
#define BAND_SOLVER 0
#define TEST_IC_CALC 0

#include <mex.h>

#include <ida/ida.h>
#include <ida/ida_direct.h> 
#define SUN_USING_SPARSE 1
#if SUNDIALS_VERSION_MAJOR >= 3
#define SUNDIALS_3 1
#endif
#if SUN_USING_SPARSE
#if ! USE_EIGEN_LU
#if SUNDIALS_3
#include <sunlinsol/sunlinsol_klu.h>
#else
#include <ida/ida_klu.h>
#endif
#endif
#elif BAND_SOLVER
#include <ida/ida_band.h>
#else
#include <ida/ida_dense.h>
#endif
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
#include "PDESolution.h"
#include "PDEEvents.h"
#include <util.h>
#include "EigenSUNSparseSolver.h"

#if SUNDIALS_3
class SunSparseMap : public FiniteDiffJacobian::SparseMap {
public:
  SunSparseMap(SUNMatrix a) :
    FiniteDiffJacobian::SparseMap(SUNSparseMatrix_Rows(a),
      SUNSparseMatrix_Columns(a), 
      SUNSparseMatrix_NNZ(a),
      SUNSparseMatrix_IndexPointers(a), 
      SUNSparseMatrix_IndexValues(a),
      SUNSparseMatrix_Data(a)) {
  }
};
#endif

PDE1dImpl::PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options) : 
pde(pde), options(options)
{
#if _cpp_lib_make_unique
  // C++14
  sf = std::make_unique<ShapeFunction2>();
#else
  sf = std::unique_ptr<ShapeFunction2>(new ShapeFunction2);
#endif
  sfm = std::unique_ptr<ShapeFunctionManager>(new ShapeFunctionManager);
  ida = 0;
  polyOrder = options.getPolyOrder();
  numIntPts = GausLegendreIntRule::getNumPtsForPolyOrder(2 * polyOrder);
  //cout << "numIntPts=" << numIntPts << endl;
  mesh = pde.getMesh();
  checkIncreasing(mesh, 5, "meshPts");
  tspan = pde.getTimeSpan();
  checkIncreasing(tspan, 6, "timePts");
  numTimes = tspan.size();
  numODE = pde.getNumODE();
  numDepVars = pde.getNumPDE();
  pdeModel = std::unique_ptr<PDEModel>(
    new PDEModel(mesh, polyOrder, numDepVars, *sfm.get()));
  numFEEqns = pdeModel->numEquations();
  totalNumEqns = numFEEqns + numODE;
  dirConsFlagsLeft.resize(numDepVars);
  dirConsFlagsRight.resize(numDepVars);

  if (numODE) {
    v.resize(numODE);
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
    isOdeAConstraint.resize(numODE);
    isOdeAConstraint.setZero();
  }
  bc.pl.resize(numDepVars);
  bc.pr.resize(numDepVars);
  bc.ql.resize(numDepVars);
  bc.qr.resize(numDepVars);

  size_t numXPts = 1;
  if (options.isVectorized()) {
    numXPts = numIntPts*pdeModel->numElements();
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
  if (options.getJacDiagnostics()) {
    cout << "Jacobian Pattern:\n";
    if (numFEEqns <= 20)
      cout << P.toDense();
    else
      cout << P;
    cout << endl;
  }
  finiteDiffJacobian = 
    std::unique_ptr<FiniteDiffJacobian>(new FiniteDiffJacobian(P));
  numViewElemsPerElem=options.getViewMesh();

  int numEvents = pde.getNumEvents();
  if (numEvents) {
    pdeEvents =
      std::unique_ptr<PDEEvents>(new PDEEvents(pde, *pdeModel));
  }
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
    cout << "resFunc called" << endl;
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
#if SUNDIALS_3
    SUNMatrix Jac, 
#else
    SlsMat Jac,
#endif
    void *user_data,
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

  int rootFunc(realtype t, N_Vector y, N_Vector yp,
    realtype *gout, void *user_data) {
    SunVector u(y);
    PDE1dImpl *fem = (PDE1dImpl*) user_data;
    fem->calcEvents(t, u, gout);
    return 0;
  }

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
  RealVector y0(totalNumEqns);
  size_t nnfe = pdeModel->numNodesFEEqns();
  getFEInitConditions(y0);
  MapMat y0FE(y0.data(), numDepVars, nnfe);
  if (numODE) {
    pde.evalODEIC(v);
    y0.bottomRows(numODE) = v;
  }
  // flag dirichlet constraints
  RealVector ul = y0FE.col(0), ur = y0FE.col(nnfe - 1);
  pde.evalBC(mesh(0), ul, mesh(mesh.size()-1), ur, tspan(0), v, vDot, bc);
  for (int i = 0; i < numDepVars; i++) {
    dirConsFlagsLeft[i] = bc.ql[i] == 0;
    dirConsFlagsRight[i] = bc.qr[i] == 0;
  }

  if (options.getEqnDiagnostics()) {
    testMats(y0);
    //sol.setSolutionVector(0, 0, y0);
    return 0;
  }
#if 0
  testODEJacobian(y0);
#endif

  const size_t neqImpl = totalNumEqns;
  /* Create vectors uu, up, res, constraints, id. */
  SunVector uu(neqImpl), up(neqImpl), res(neqImpl), id(neqImpl);
  double *udata = &uu[0];
  MapVec u(udata, neqImpl);
  std::copy_n(y0.data(), totalNumEqns, udata);
  up.setConstant(0);
#if 0
  // test for consistent initial conditions
  resFunc(0, uu(), up(), res(), this);
  MapVec resVec(&res[0], neqImpl);
  double initResNorm = resVec.dot(resVec);
  pdePrintf("initResNorm=%12.3e\n", sqrt(initResNorm));
#endif

  /* Call IDACreate and IDAMalloc to initialize solution */
  ida = IDACreate();
  check_flag(ida, "IDACreate", 0);

  int ier = IDASetUserData(ida, this);
  check_flag(&ier, "IDASetUserData", 1);
  if (!options.getICMethod() || numODE) {
    setAlgVarFlags(uu, up, id);
    ier = IDASetId(ida, id.getNV());
    check_flag(&ier, "IDASetId", 1);
  }
  if (numODE) {
    // IDASetId needed for this call
    ier = IDASetSuppressAlg(ida, true);
    check_flag(&ier, "IDASetSuppressAlg", 1);
  }
  double t0 = tspan(0);
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
#if SUNDIALS_3
  SUNMatrix A = SUNSparseMatrix((sunindextype)totalNumEqns,
    (sunindextype)totalNumEqns, (sunindextype)numNonZerosJacMax, CSC_MAT);
  check_flag(A, "SUNSparseMatrix", 0);
#if USE_EIGEN_LU
  //cout << "Using Eigen Sparse LU" << endl;
  EigenSUNSparseSolver linearSolver(uu.getNV(), A);
  SUNLinearSolver LS = &linearSolver;
#else
  SUNLinearSolver LS = SUNKLU(uu.getNV(), A);
  check_flag(LS, "SUNKLU", 0);
#endif
  ier = IDADlsSetLinearSolver(ida, LS, A);
  check_flag(&ier, "IDADlsSetLinearSolver", 1);
  ier = IDADlsSetJacFn(ida, jacFunc);
  check_flag(&ier, "IDADlsSetJacFn", 1);
#else
  ier = IDAKLU(ida, (int) totalNumEqns, (int) numNonZerosJacMax, CSC_MAT);
  check_flag(&ier, "IDAKLU", 1);
  ier = IDASlsSetSparseJacFn(ida, jacFunc);
  check_flag(&ier, "IDASlsSetSparseJacFn", 1);
#endif
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
  uu = initCond.getU0();
  up = initCond.getUp0();

#if TEST_IC_CALC
  // testing only
  SunVector u0Tmp = initCond.getU0();
  SunVector up0Tmp = initCond.getUp0();
  testICCalc(u0Tmp, up0Tmp, res, id, tf);
#endif

  int numEvents = pde.getNumEvents();
  if (numEvents) {
    ier = IDARootInit(ida, numEvents, rootFunc);
    check_flag(&ier, "IDARootInit", 1);
  }
  
  //sol.time(0) = tspan(0);
  sol.setSolutionVector(0, 0, initCond.getU0().topRows(numFEEqns));
  if (numODE) {
    sol.uOde.resize(numTimes, numODE);
    sol.uOde.row(0) = initCond.getU0().bottomRows(numODE);
  }

  // optionally, calc and print jacobian matrices
  if(options.getJacDiagnostics())
     jacobianDiagnostics(tspan(0), uu, up, res);

  bool doTerm = false;
  int i = 1;
  while (i< numTimes && ! doTerm) {
    double tout=tspan(i), tret;
    ier = IDASolve(ida, tout, &tret, uu.getNV(), up.getNV(), IDA_NORMAL);
    if (ier < 0) {
      pdePrintf("Error returned from IDASolve=%d\n", ier);
      printStats();
      char msg[1024];
      sprintf(msg, "Time integration failed at t=%15.6e before reaching final time.\n"
        "Often this is caused by one or more dependent variables becoming unbounded.",
        tret);
      throw PDE1dException("pde1d:integ_failure", msg);
    }
    //cout << "tret=" << tret << endl;
    // events
    bool eventSatisfied = false;
    if (ier == IDA_ROOT_RETURN) {
      eventSatisfied = true;
      IntVector eventsFound(numEvents);
      ier = IDAGetRootInfo(ida, eventsFound.data());
      check_flag(&ier, "IDAGetRootInfo", 1);
      doTerm = pdeEvents->isTerminalEvent(tret, u, eventsFound);
      sol.setEventsSolution(i, tret, u.topRows(numFEEqns), eventsFound);
    }
    if (!eventSatisfied || doTerm) {
      sol.setSolutionVector(i, tret, u.topRows(numFEEqns));
      if (numODE)
        sol.uOde.row(i) = u.bottomRows(numODE);
      ++i;
    }
  }
  printStats();
  sol.close();

  return 0;
}

void PDE1dImpl::calcGlobalEqns(double time, SunVector &u, SunVector &up,
  RealVector &Cxd, RealVector &F, RealVector &S)
{
  Cxd.setZero(); F.setZero(); S.setZero();
  if (options.isVectorized() && pde.hasVectorPDEEval())
    calcGlobalEqnsVectorized(time, u, up, Cxd, F, S);
  else
    calcGlobalEqnsNonVectorized(time, u, up, Cxd, F, S);
}

template<class T, class TR>
void PDE1dImpl::calcGlobalEqnsNonVectorized(double t, T &u, T &up,
  TR &Cxd, TR &F, TR &S)
{
  bool useDiagMassMat = options.getDiagMassMat();
  const ShapeFunctionManager::EvaluatedSF &esf =
    sfm->getShapeFunction(polyOrder);
  const RealMatrix &N = esf.N();
  const RealMatrix &dN = esf.dN();
  const RealVector &intWts = esf.intRuleWts();

  const size_t nnfee = pdeModel->numNodesFEEqns();
  const size_t nen = N.rows();
  const int m = pde.getCoordSystem();
  size_t numElemEqns = numDepVars*nen;
  Cxd.setZero();

  //cout << "u=" << u.transpose() << endl;
  MapMat u2(u.data(), numDepVars, nnfee);
  MapMat up2(up.data(), numDepVars, nnfee);
  //cout << "u2=" << u2.transpose() << endl;
  //cout << "up2=" << up2.transpose() << endl;
  RealMatrix u2e(numDepVars, nen);

#if 0
  cout << "intWts=" << intWts.transpose() << endl;
  cout << "N\n" << N << endl;
#endif

  pdeCoeffs.c.resize(numDepVars, 1);
  pdeCoeffs.f.resize(numDepVars, 1);
  pdeCoeffs.s.resize(numDepVars, 1);

  RealVector eF = RealVector::Zero(numElemEqns);
  MapMat eFX(eF.data(), numDepVars, nen);
  RealVector eS = RealVector::Zero(numElemEqns);
  MapMat eSX(eS.data(), numDepVars, nen);
  RealVector eC = RealVector::Zero(numElemEqns);
  MapMat eCX(eC.data(), numDepVars, nen);
  RealMatrix up2e(numDepVars, nen);
  Eigen::Map<RealVector> eUp(up2e.data(), numElemEqns);
  PDEModel::DofList eDofs(nen);
  RealMatrix oNen = RealMatrix::Ones(1,nen);
  RealVector dNdx(nen), ui(numDepVars), upi(numDepVars), dUiDx(numDepVars);

  RealVector &x = mesh;
  size_t ne = pdeModel->numElements();
  for (int e = 0; e < ne; e++) {
    pdeModel->getDofIndicesForElem(e, eDofs);
    PDEModel::globalToElemVec(eDofs, u2, u2e);
    PDEModel::globalToElemVec(eDofs, up2, up2e);
    eC.setZero();
    eF.setZero();
    eS.setZero();
    double L = x(e + 1) - x(e);
    double jac = L / 2;
    for (int i = 0; i < numIntPts; i++) {
      // assume two node elements
      double jacWt = jac*intWts(i);
      double xi = x(e)*N(0, i) + x(e + 1)*N(1, i);
      dNdx = dN.col(i) / jac;
      ui = u2e*N.col(i);
      upi = up2e*N.col(i);
      dUiDx = u2e*dNdx;
      pde.evalPDE(xi, t, ui, dUiDx, v, vDot, pdeCoeffs);
      bool isFullCMat = numDepVars>1 && pdeCoeffs.c.rows() == pdeCoeffs.c.cols();
      Eigen::Ref<Eigen::VectorXd> sIp = pdeCoeffs.s.col(0);
      Eigen::Ref<Eigen::VectorXd> fIp = pdeCoeffs.f.col(0);
      checkCoeffs(pdeCoeffs);
      double xm = 1;
      if (m == 1)
        xm = xi;
      else if(m == 2)
        xm = xi*xi;

      pdeCoeffs.c *= xm;
      sIp *= xm;
      fIp *= xm;
 
      eSX += sIp * N.col(i).transpose() * jacWt;
      eFX += fIp * dNdx.transpose() * jacWt;
      if(isFullCMat)
        eCX += pdeCoeffs.c* upi* N.col(i).transpose() * jacWt;
      else {
        if (useDiagMassMat) {
          eCX += pdeCoeffs.c * oNen * jacWt / (double)nen; // lump mass
        }
        else
          eCX += (pdeCoeffs.c.array() * upi.array()).matrix() * N.col(i).transpose() * jacWt;
      }
#if 0
      printMat(upi, "upi");
      printMat(N.col(i).transpose(), "N");
      printMat(eCX, "eCX");
#endif
    } // end integration point loop
#if DEBUG_MATS
    cout << "eF: " << eF.transpose() << endl;
    if (useDiagMassMat)
      cout << "eCMatD=" << eCMatD.diagonal().transpose() << endl;
    else
      cout << "eCMat\n" << eCMat << endl;
#endif
    if (useDiagMassMat) {
      eCX = eCX.array() * up2e.array();
    }
    PDEModel::assembleElemVec(eDofs, nnfee, numDepVars, eC, Cxd);
    //printMat(eC, "eC");
    //printMat(eCX, "eCX");
    PDEModel::assembleElemVec(eDofs, nnfee, numDepVars, eF, F);
    PDEModel::assembleElemVec(eDofs, nnfee, numDepVars, eS, S);
  }
}

  template<class T, class TR>
  void PDE1dImpl::calcGlobalEqnsVectorized(double t, T &u, T &up,
    TR &Cxd, TR &F, TR &S)
  {
    bool useDiagMassMat = options.getDiagMassMat();
    const ShapeFunctionManager::EvaluatedSF &esf =
      sfm->getShapeFunction(polyOrder);
    const RealMatrix &N = esf.N();
    const RealMatrix &dN = esf.dN();
    const RealVector &intWts = esf.intRuleWts();

    const size_t nnfee = pdeModel->numNodesFEEqns();
    const size_t nen = N.rows();
    const int m = pde.getCoordSystem();
    size_t numElemEqns = numDepVars*nen;
    Cxd.setZero();

    //cout << "u=" << u.transpose() << endl;
    MapMat u2(u.data(), numDepVars, nnfee);
    MapMat up2(up.data(), numDepVars, nnfee);
    //cout << "u2=" << u2.transpose() << endl;
    //cout << "up2=" << up2.transpose() << endl;
    RealMatrix u2e(numDepVars, nen);
    RealVector dNdx(nen), ui(numDepVars), upi(numDepVars), dUiDx(numDepVars);

#if 0
    cout << "intWts=" << intWts.transpose() << endl;
    cout << "N\n" << N << endl;
#endif

    // get the coefficients at all integ pts in a single call
    size_t ne = pdeModel->numElements();
    size_t numXPts = numIntPts*ne;
    xPts.resize(numXPts);
    uPts.resize(numDepVars, numXPts);
    duPts.resize(numDepVars, numXPts);

    PDEModel::DofList eDofs(nen);

    RealVector &x = mesh;
    int ip = 0;
    for (int e = 0; e < ne; e++) {
      pdeModel->getDofIndicesForElem(e, eDofs);
      PDEModel::globalToElemVec(eDofs, u2, u2e);
      //cout << "u2e=" << u2e << endl;
      double L = x(e + 1) - x(e);
      double jac = L / 2;
      for (int i = 0; i < numIntPts; i++) {
        xPts(ip) = x(e)*N(0, i) + x(e + 1)*N(1, i);
        dNdx = dN.col(i) / jac;
        uPts.col(ip) = u2e*N.col(i);
        duPts.col(ip) = u2e*dNdx;
        ip++;
      }
    }
    pdeCoeffs.c.resize(numDepVars, numXPts);
    pdeCoeffs.f.resize(numDepVars, numXPts);
    pdeCoeffs.s.resize(numDepVars, numXPts);
    pde.evalPDE(xPts, t, uPts, duPts, v, vDot, pdeCoeffs);

    RealVector eF = RealVector::Zero(numElemEqns);
    MapMat eFX(eF.data(), numDepVars, nen);
    RealVector eS = RealVector::Zero(numElemEqns);
    MapMat eSX(eS.data(), numDepVars, nen);
    RealVector eC = RealVector::Zero(numElemEqns);
    MapMat eCX(eC.data(), numDepVars, nen);
    RealMatrix up2e(numDepVars, nen);
    Eigen::Map<RealVector> eUp(up2e.data(), numElemEqns);
    RealMatrix oNen = RealMatrix::Ones(1, nen);

    ip = 0; 
    for (int e = 0; e < ne; e++) {
      pdeModel->getDofIndicesForElem(e, eDofs);
      PDEModel::globalToElemVec(eDofs, up2, up2e);
      eF.setZero();
      eS.setZero();
      eC.setZero();
      double L = x(e + 1) - x(e);
      double jac = L / 2;
      for (int i = 0; i < numIntPts; i++) {
        // assume two node elements
        double jacWt = jac*intWts(i);
        double xi = xPts(ip);
        dNdx = dN.col(i) / jac;
        upi = up2e * N.col(i);
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

        eSX += sIp * N.col(i).transpose() * jacWt;
        eFX += fIp * dNdx.transpose() * jacWt;
        if (useDiagMassMat)
          eCX += cIp * oNen * jacWt / (double)nen; // lump mass
        else
          eCX += (cIp.array() * upi.array()).matrix() * N.col(i).transpose() * jacWt;
#if 0
        printMat(cIp, "cIp");
        printMat(upi, "upi");
        printMat(N.col(i).transpose(), "N");
        printMat(eCX, "eCX");
#endif
        ip++;
      } // end integration point loop
#if DEBUG_MATS
      cout << "eF: " << eF.transpose() << endl;
      cout << "eCMat\n" << eCMat << endl;
#endif
      if (useDiagMassMat) {
        eCX = eCX.array() * up2e.array();
      }
      PDEModel::assembleElemVec(eDofs, nnfee, numDepVars, eC, Cxd);
      //cout << "eF=" << eF.transpose() << endl;
      //printMat(eC, "eC");
      PDEModel::assembleElemVec(eDofs, nnfee, numDepVars, eF, F);
      PDEModel::assembleElemVec(eDofs, nnfee, numDepVars, eS, S);
    }
#if 0
    cout << "F\n" << F.transpose() << endl;
    cout << "C\n" << Cxd.transpose() << endl;
    cout << "S\n" << S.transpose() << endl;
#endif
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

    const size_t nnfe = pdeModel->numNodesFEEqns();
    MapMat u2(u.data(), numDepVars, nnfe);

    // add odes, if any
    if (numODE) {
      MapMat f2(F.data(), numDepVars, nnfe);
      meshMapper->mapFunction(u2, odeU);
      meshMapper->mapFunctionDer(u2, odeDuDx);
      meshMapper->mapFunction(f2, odeFlux);
      MapMat up2(up.data(), numDepVars, nnfe);
      meshMapper->mapFunction(up2, odeDuDt);
      meshMapper->mapFunctionDer(up2, odeDuDxDt);
      pde.evalODE(time, v, vDot, odeU, odeDuDx, odeFlux,
        odeDuDt, odeDuDxDt, odeF);
      R.bottomRows(numODE) = odeF;

      // add lagrange multiplier terms
#if 0
      RealMatrix dOduDu = calcDOdeDu(time, odeU,
        ypFE, r2, v, vdot);
#endif
    }

  // apply constraints
  size_t rightDofOff = numFEEqns - numDepVars;
  RealVector ul = u2.col(0), ur = u2.col(nnfe-1);
  const double xl = mesh(0), xr = mesh(mesh.size() - 1);
  pde.evalBC(xl, ul, xr, ur, time, v, vDot, bc);
  const int m = pde.getCoordSystem();
  bool sing = m > 0 && xl == 0;
  for (int i = 0; i < numDepVars; i++) {
    if (bc.ql(i) != 0) {
      //printf("left bc: i=%d, ql=%f, pl=%f\n", i, bc.ql(i), bc.pl(i));
      double qli = bc.ql(i);
      if (!sing) {
        if (m == 1)
          qli /= xl;
        else if (m == 2)
          qli /= (xl * xl);
      }
      R(i) -= bc.pl(i) / qli;
    }
    else {
      // apply dirichlet constraint
      R(i) = -bc.pl(i);
      Cxd(i) = 0;
    }
    if (bc.qr(i) != 0) {
      double qri = bc.qr(i);
      if (m == 1)
        qri /= xr;
      else if (m == 2)
        qri /= (xr*xr);
      R(i + rightDofOff) += bc.pr(i) / qri;
    }
    else {
      // apply dirichlet constraint
      R(i + rightDofOff) = -bc.pr(i);
      Cxd(i + rightDofOff) = 0;
    }
  }
  R.topRows(numFEEqns) += Cxd;
#if 0
  cout << "Cxd\n" << Cxd << endl;
  cout << "F\n" << F << endl;
#endif
  
}

  void PDE1dImpl::jacobianDiagnostics(double t0, SunVector &u,
     SunVector &up, SunVector &R)
  {
    int jacDiag = options.getJacDiagnostics();
    if (!jacDiag)
      return;

    R.setZero();
    calcRHSODE(0, u, up, R);
    cout << "Residual:\n" << R.transpose() << endl;

    bool printDense = numFEEqns <= 12;
    SparseMat jac(numFEEqns, numFEEqns);
    calcJacobian(t0, 1, 0, u, up, R, jac);
    cout << "DfDu:\n";
    if(printDense)
      cout << jac.toDense();
    else
      cout << jac;
    cout << endl;

    calcJacobian(t0, 0, 1, u, up, R, jac);
    cout << "DfDuDot:\n";
    if (printDense)
      cout << jac.toDense();
    else
      cout << jac;
    cout << endl;

    if(abs(jacDiag) > 1)
    throw PDE1dException("pde1d:testJacobian", "Test complete.");
  }

  PDE1dImpl::MatrixVec PDE1dImpl::testODEJacobian(RealVector &y) {
    MatrixVec mats;
    if (numODE) {
      double t0 = tspan[0];
      SunVector yTmp(totalNumEqns), yp(totalNumEqns);
      yTmp = y;
      yp.setZero();
      calcGlobalEqns(t0, yTmp, yp, Cxd, F, S);
      const size_t nnfe = pdeModel->numNodesFEEqns();
      MapMat  f2(F.data(), numDepVars, nnfe);
      MapMat  yFE(y.data(), numDepVars, nnfe);
      MapMat  ypFE(yp.data(), numDepVars, nnfe);
      v = y.bottomRows(numODE);
      vDot = yp.bottomRows(numODE);
      RealMatrix dOdeDvDot = calcDOdeDvDot(t0, yFE, ypFE, f2, v, vDot);
      RealMatrix dOdeDv = calcDOdeDv(t0, yFE, ypFE, f2, v, vDot);
      cout << "dOdeDv:\n" << dOdeDv << endl;
      cout << "dOdeDvDot:\n" << dOdeDvDot << endl;

      RealMatrix dOdeDu, dOdeDuDot;
      calcDOdeDu(t0, yFE, ypFE, f2, v, vDot, dOdeDu, dOdeDuDot);
      cout << "dOdeDu:\n" << dOdeDu << endl;
      cout << "dOdeDuDot:\n" << dOdeDuDot << endl;
      mats = { dOdeDv, dOdeDvDot, dOdeDu, dOdeDuDot };
    }
    return mats;
    //throw PDE1dException("pde1d:testODEJacobian", "Test complete.");
  }

void PDE1dImpl::testMats(const RealVector &y0)
{
  const int diag = options.getEqnDiagnostics();
  pdePrintf("Begin testMats: numFEEqns = %d\n", numFEEqns);
  cout << "totalNumEqns=" << totalNumEqns << endl;
  SunVector u(totalNumEqns), up(totalNumEqns), R(totalNumEqns);
#if 0
  for (int i = 0; i < totalNumEqns; i++) {
    u(i) = i;
    up(i) = i;
  }
#elif 1
  u = y0;
  double upEnd = static_cast<double>(totalNumEqns - 1);
  up = Eigen::VectorXd::LinSpaced(totalNumEqns, 0, upEnd);
#else
  u = y0;
  up.setOnes();
#endif
#if 1
  printSystemVector(u, "u");
  printSystemVector(up, "up");
  RealVector Cxd = RealVector::Zero(totalNumEqns);
  RealVector F = RealVector::Zero(totalNumEqns);
  RealVector S = RealVector::Zero(totalNumEqns);
  // copy the ode dofs to their own vectors
  if (numODE) {
    v = u.bottomRows(numODE);
    vDot = up.bottomRows(numODE);
  }
  calcGlobalEqns(0, u, up, Cxd, F, S);
  printSystemVector(Cxd, "Cxd");
  printSystemVector(F, "F");
  printSystemVector(S, "S");
#if 1
  R.setZero();
  calcRHSODE(0, u, up, R);
#else
  R = S - F;
#endif
  printSystemVector(R, "R");
#endif
  if (diag > 1) {
    SparseMat jac(totalNumEqns, totalNumEqns);
    calcJacobian(0, 1, 0, u, up, R, jac);
    cout << "jac\n" << jac.toDense() << endl;
  }
  if (diag > 1) {
    SparseMat jac(totalNumEqns, totalNumEqns);
    calcJacobian(0, 0, 1, u, up, R, jac);
    cout << "mass matrix\n" << jac.toDense() << endl;
    if (options.getDiagMassMat()) {
      Eigen::VectorXd mDiag = jac.diagonal();
      printSystemVector(mDiag, "M-diagonal");
    }
  }
  throw PDE1dException("pde1d:testResidual", "test complete.");
}


void PDE1dImpl::setAlgVarFlags(SunVector &y0, SunVector &y0p, SunVector &id)
{
  const size_t nnfe = pdeModel->numNodesFEEqns();
  MapMat y0FE(y0.data(), numDepVars, nnfe);
  
  double t0 = tspan[0];

  id.setConstant(1);
  double *iddata = id.data();
  MapMat idMat(iddata, numDepVars, nnfe);

  // flag equations where c==0
  double xi[] = { -1, 1 };
  const int nen = 2;
  RealMatrix dN(nen, nen);
  for (int i = 0; i < nen; i++) {
    double pt = xi[i];
    sf->dNdr(pt, dN.col(i).data());
  }
 
  const size_t nnMesh = mesh.size();
  duPts.resize(numDepVars, nnMesh);
  RealMatrix uPts(numDepVars, nnMesh);
  size_t nnm1 = nnMesh - 1;
  size_t iu = 0;
  for (int i = 0; i < nnMesh; i++) {
    double L;
    if (i < nnm1) {
      L = mesh(i + 1) - mesh(i);
    }
    else
      L = mesh(i) - mesh(i - 1);
    double jac = L / 2;
    auto dNdx = dN.col(0) / jac;
    uPts.col(i) = y0FE.col(iu);
    if (i < nnm1) {
      duPts.col(i) = (y0FE.col(iu)*dNdx(0) + y0FE.col(iu + 1)*dNdx(1));
      int  nen = pdeModel->element(i).numNodes();
      iu += nen - 1;
    }
    else
      duPts.col(i) = (y0FE.col(iu-1)*dNdx(0) + y0FE.col(iu)*dNdx(1));
  }

  size_t numXPts = 1;
  if (options.isVectorized())
    numXPts = nnMesh;
  pdeCoeffs.c.resize(numDepVars, numXPts);
  pdeCoeffs.f.resize(numDepVars, numXPts);
  pdeCoeffs.s.resize(numDepVars, numXPts);

  if (options.isVectorized() && pde.hasVectorPDEEval()) {
    pde.evalPDE(mesh, 0, uPts, duPts, 
      v, vDot, pdeCoeffs);
    for (int i = 0; i < nnMesh; i++) {
      for (int j = 0; j < numDepVars; j++)
        if (pdeCoeffs.c(j,i)==0)
          idMat(j, i) = 0;
    }
  }
  else {
    for (int i = 0; i < nnMesh; i++) {
      double xi = mesh(i);
      const auto &ui = uPts.col(i);
      const auto &dUiDx = duPts.col(i);
      Eigen::Ref<Eigen::VectorXd> cIp = pdeCoeffs.c.col(0);
      pde.evalPDE(xi, t0, ui, dUiDx, v, vDot, pdeCoeffs);
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
    const size_t nnfe = pdeModel->numNodesFEEqns();
    MapMat  f2(F.data(), numDepVars, nnfe);
    MapMat  yp0FE(y0p.data(), numDepVars, nnfe);
    v = y0.bottomRows(numODE);
    vDot = y0p.bottomRows(numODE);
    RealMatrix odeJacDot = calcDOdeDvDot(t0, y0FE, yp0FE, f2, v, vDot);
    //cout << "odeJac:\n" << odeJac << endl;
    for (int i = 0; i < numODE; i++) {
      auto rci = odeJacDot.row(i) + odeJacDot.col(i).transpose();
      if ((rci.array() == 0).all())
        id(numFEEqns + i) = 0;
    }
    //print(id.getNV(), "id");

    // Do the ODE equations depend explicitly on their respective ODE
    // variable? If not, the equation is a pure constraint on the values
    // of the PDE variables
    RealMatrix odeJac = calcDOdeDv(t0, y0FE, yp0FE, f2, v, vDot);
    for (int i = 0; i < numODE; i++) {
      //pdePrintf("odeJac.all=%d\n", (odeJac.row(i).array() == 0).all());
      //pdePrintf("odeJacDot.all=%d\n", (odeJacDot.row(i).array() == 0).all());
      if ((odeJacDot.row(i).array() == 0).all() &&
        (odeJac.row(i).array() == 0).all())
        isOdeAConstraint(i) = 1;
      //pdePrintf("isOdeAConstraint=%d\n", i);
    }
    //printMat(odeJac, "odeJac", "%12.3e");
    //printMat(odeJacDot, "odeJacDot", "%12.3e");
    //printMat(isOdeAConstraint.transpose(), "####isOdeAConstraint", "%2d");
  }
}

RealMatrix PDE1dImpl::calcDOdeDvDot(double time, const RealMatrix &yFE,
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

RealMatrix PDE1dImpl::calcDOdeDv(double time, const RealMatrix &yFE,
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
  RealVector vDelta = v, fDelta(numODE);
  double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
  for (int i = 0; i < numODE; i++) {
    double h = sqrtEps * std::max(v[i], 1.0);
    double tmpVi = v[i];
    vDelta[i] += h;
    pde.evalODE(time, vDelta, vdot, odeU, odeDuDx, odeFlux,
      odeDuDt, odeDuDxDt, fDelta);
    jac.col(i) = (fDelta - odeF) / h;
    vDelta[i] = tmpVi;
  }
  return jac;
}

void PDE1dImpl::calcDOdeDu(double time, const RealMatrix &yFE,
  const RealMatrix &ypFE, const RealMatrix &r2, RealVector &v,
  RealVector &vdot, RealMatrix &jac, RealMatrix &jacDot)
{
  // calculate derivatives of ode equations wrt pde variables
  const PDEMeshMapper::StdIntVec &pdeDofs = meshMapper->mappedDOFList();
  const size_t numPdeDofs = pdeDofs.size();
  if (!numPdeDofs) return;
  cout << "pdeDofs=" << pdeDofs << endl;
  jac.resize(numODE, numPdeDofs*numDepVars);
  jacDot.resize(numODE, numPdeDofs*numDepVars);

  meshMapper->mapFunction(yFE, odeU);
  meshMapper->mapFunctionDer(yFE, odeDuDx);
  meshMapper->mapFunctionDer(r2, odeFlux);
  meshMapper->mapFunction(ypFE, odeDuDt);
  meshMapper->mapFunctionDer(ypFE, odeDuDxDt);
  pde.evalODE(time, v, vDot, odeU, odeDuDx, odeFlux,
    odeDuDt, odeDuDxDt, odeF);

  double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
  // dODE/dU
  RealMatrix uDelta = yFE;
  RealVector fDelta(numODE);
  int k = 0;
  for (int i = 0; i < numPdeDofs; i++) {
    int ui = pdeDofs[i];
    for (int j = 0; j < numDepVars; j++) {
      double tmpUij = yFE(j, ui);
      double h = sqrtEps * std::max(tmpUij, 1.0);
      uDelta(j, ui) += h;
      meshMapper->mapFunction(uDelta, odeU);
      meshMapper->mapFunctionDer(uDelta, odeDuDx);
      pde.evalODE(time, v, vdot, odeU, odeDuDx, odeFlux,
        odeDuDt, odeDuDxDt, fDelta);
      jac.col(k++) = (fDelta - odeF) / h;
      uDelta(j, ui) = tmpUij;
    }
  }

  // dODE/dUDot
  meshMapper->mapFunction(yFE, odeU);
  meshMapper->mapFunctionDer(yFE, odeDuDx);

  RealMatrix dudtDelta = ypFE;
  k = 0;
  for (int i = 0; i < numPdeDofs; i++) {
    int ui = pdeDofs[i];
    for (int j = 0; j < numDepVars; j++) {
      double tmpUij = ypFE(j, ui);
      double h = sqrtEps * std::max(tmpUij, 1.0);
      dudtDelta(j, ui) += h;
      meshMapper->mapFunction(dudtDelta, odeDuDt);
      meshMapper->mapFunctionDer(dudtDelta, odeDuDxDt);
      pde.evalODE(time, v, vdot, odeU, odeDuDx, odeFlux,
        odeDuDt, odeDuDxDt, fDelta);
      jacDot.col(k++) = (fDelta - odeF) / h;
      dudtDelta(j, ui) = tmpUij;
    }
  }
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
  pdePrintf("Number of internal time steps = %ld\n", nsteps);
  pdePrintf("Number of residual function calls = %ld\n", nrevals);
  pdePrintf("Number of Jacobian calculations = %ld\n", nlinsetups);
  pdePrintf("Number of solution accuracy test failures = %ld\n", netfails);
  pdePrintf("Last internal time step size = %12.3e\n", hlast);
  pdePrintf("Number of nonlinear iterations = %ld\n", nniters);
  pdePrintf("Number of nonlinear convergence failures = %ld\n", nncfails);
}

void PDE1dImpl::calcJacPattern(Eigen::SparseMatrix<double> &J)
{
  J.resize(totalNumEqns, totalNumEqns);
  size_t n2 = numDepVars*numDepVars;
  size_t nel = pdeModel->numElements();
#if 1
  size_t nnz = 0;
  size_t maxNN = 0;
  for (int i = 0; i < nel; i++) {
    size_t nen = pdeModel->element(i).numNodes();
    maxNN = std::max(maxNN, nen);
    nnz += nen*nen*n2;
  }
  typedef Eigen::Triplet<int, int> T;
  std::vector<T> tripList;
  tripList.reserve(nnz);
  int eOff = 0;
  for (int i = 0; i < nel; i++) {
    int nen = pdeModel->element(i).numNodes();
    int neleq = nen*(int)numDepVars;
    for (int i = 0; i < neleq; i++) {
      for (int j = 0; j < neleq; j++) {
        tripList.push_back(T(i+eOff, j+eOff, 1));
      }
    }
    eOff += (nen-1)*static_cast<int>(numDepVars);
  }
  J.setFromTriplets(tripList.begin(), tripList.end());
#else
  size_t nnz = 3 * n2*(nel + 1); // approximate nnz
  //cout << "nnz=" << nnz << ", nnzX=" << nnzX << endl;
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
#endif
  // odes may couple with any other vars
  eOff = static_cast<int>(numFEEqns);
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

#if SUNDIALS_3
void PDE1dImpl::calcJacobianODE(double time, double beta, SunVector &u, 
  SunVector &up, SunVector &res, SUNMatrix Jac)
{
  SunSparseMap eigJac(Jac);
#else
void PDE1dImpl::calcJacobianODE(double time, double beta, SunVector &u,
  SunVector &up, SunVector &res, SlsMat Jac)
{
FiniteDiffJacobian::SparseMap eigJac(Jac->N, Jac->M,
  Jac->NNZ, Jac->indexptrs, Jac->indexvals, Jac->data);
#endif
  const bool useCD = !true; // use central difference approximation, if true
  finiteDiffJacobian->calcJacobian(time, 1, beta, u.getNV(), up.getNV(),
    res.getNV(), resFunc, this, eigJac, useCD);
#if 0
#if SUNDIALS_3
  SUNSparseMatrix_Print(Jac, stdout);
#else
  SparsePrintMat(Jac, stdout);
#endif
#endif
}

void PDE1dImpl::calcJacobian(double time, double alpha, double beta, SunVector &u,
  SunVector &up, SunVector &R, SparseMat &Jac)
{
  const bool useCD = !true; // use central difference approximation, if true
  finiteDiffJacobian->calcJacobian(time, alpha, beta, u.getNV(), up.getNV(), 
    R.getNV(), resFunc, this, Jac, useCD);
}

void PDE1dImpl::calcEvents(double time, const SunVector &u, 
  double *gOut)
{
  pdeEvents->calcEvents(time, u, gOut);
}

#if TEST_IC_CALC
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
  printf("dfdy: n=%td, rank=%td\n", dfDy.rows(), lu.rank());

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
#endif

double PDE1dImpl::calcResidualNorm(double t, SunVector &uu, SunVector &up, SunVector &res)
{
  resFunc(t, uu.getNV(), up.getNV(), res.getNV(), this);
  MapVec r(&res[0], res.size());
  double rRms = sqrt(r.dot(r));
  return rRms;
}

void PDE1dImpl::getFEInitConditions(RealVector &y0)
{
  size_t nnfe = pdeModel->numNodesFEEqns();
  MapMat y0FE(y0.data(), numDepVars, nnfe);
  RealVector ic(numDepVars);
  const int md = ShapeFunctionHierarchical::MAX_DEGREE;
  const int mdp1 = md + 1;
  RealVector xm(mdp1);
  RealMatrix Nx(mdp1, mdp1), Ux(mdp1, numDepVars), Ui(mdp1, numDepVars);
  PDEModel::DofList dofs(mdp1);
  Eigen::PartialPivLU<RealMatrix> lu;
  int nnLast = 0;
  int iOff = 0;
  for (int i = 0; i < pdeModel->numElements(); i++) {
    const PDEElement &elem = pdeModel->element(i);
    int nn = elem.numNodes();
    xm = RealVector::LinSpaced(nn, mesh(i), mesh(i+1));
    int start = 1;
    if (!i) start = 0;
    if (nn == 2) {
      for (int j = start; j < nn; j++) {
        pde.evalIC(xm(j), ic);
        y0FE.col(iOff++) = ic;
      }
    }
    else {   
      pdeModel->getDofIndicesForElem(i, dofs);  
      Ux.resize(nn, numDepVars);
      for (int j = 0; j < nn; j++) {
        pde.evalIC(xm(j), ic);
        Ux.row(j) = ic;
      }
      //cout << "Ux\n" << Ux << endl;
      if (nn != nnLast) {
        Nx.resize(nn, nn);
        const ShapeFunction &sf = elem.getSF().getShapeFunction();
        double r = -1, dr = 2. / (double)(nn - 1);
        for (int j = 0; j < nn; j++) {
          sf.N(r, Nx.col(j).data());
          r += dr;
        }
        //cout << "Nx\n" << Nx << endl;
        lu.compute(Nx.transpose());
        nnLast = nn;
      }
      Ui = lu.solve(Ux);
      //cout << "Ui\n" << Ui << endl;
      PDEModel::elemToGlobalVec(dofs, Ui.transpose(), y0FE);
    }
  }
  //cout << "y0FE\n" << y0FE << endl;
}

template<class T>
void PDE1dImpl::interpolateGlobalVecToViewMesh(const T &gVec,
  RealMatrix &uViewNodes)
{
  const size_t nnfee = pdeModel->numNodesFEEqns();
  typedef Eigen::Map<const Eigen::MatrixXd> ConstMapMat;
  ConstMapMat u2(gVec.data(), numDepVars, nnfee);
  //RealVector &x = mesh;
  PDEModel::DofList eDofs;
  RealMatrix u2e(numDepVars, 2);
  size_t numR = numViewElemsPerElem - 1;
  double dr = 2. / (double)numViewElemsPerElem;
  RealVector rVals(numR);
  double r = -1;
  for (int i = 0; i < numR; i++) {
    r += dr;
    rVals(i) = r;
  }
  cout << "rvals=" << rVals << endl;
  RealMatrix N, u2eView;
  size_t ne = pdeModel->numElements();
  int lastNumElemNodes = 0;
  size_t iViewNode = 0;
  for (int e = 0; e < ne; e++) {
    const PDEElement &elem = pdeModel->element(e);
    pdeModel->getDofIndicesForElem(e, eDofs);
    PDEModel::globalToElemVec(eDofs, u2, u2e);
    cout << "u2e=" << u2e << endl;
    const ShapeFunction &sf = elem.getSF().getShapeFunction();
    if (sf.numNodes() != lastNumElemNodes) {
      lastNumElemNodes = sf.numNodes();
      N.resize(lastNumElemNodes, numR);
      for(int i=0; i<numR; i++)
        sf.N(rVals(i), N.col(i).data());
      cout << "N\n" << N << endl;
    }
    u2eView = u2e*N;
    //cout << "uViewNodes\n" << uViewNodes << endl;
    uViewNodes.col(iViewNode) = u2e.col(0);
#if 0
    for (int i = 0; i < numR; i++) {
      uViewNodes.col(iViewNode + i + 1) = u2eView.col(i);
    }
#else
    uViewNodes.block(0, iViewNode + 1, numDepVars, numR) = u2eView;
#endif
    iViewNode += numR + 1;
  }
  // finish the last view node
  uViewNodes.col(iViewNode) = u2e.col(1);
}

template<class T>
void PDE1dImpl::printSystemVector(const T& v, const char* name) {
  //cout << name << '\n';
  if (!numODE) {
    auto np = numFEEqns / numDepVars;
    Eigen::Map<const Eigen::MatrixXd> mv(v.data(), numDepVars, np);
    cout << mv << '\n';
    printMat(mv, name);
  }
  else
    //cout << v.transpose() << '\n';
    printMat(v.transpose(), name);
}
