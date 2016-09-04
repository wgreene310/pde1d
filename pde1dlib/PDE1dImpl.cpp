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

using std::cout;
using std::endl;

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
#if SUN_USING_SPARSE
#include "FiniteDiffJacobian.h"
#endif

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

#define DEBUG_MATS 0

namespace {
  // shape functions and derivatives for 2-node element
  void shapeLine2(double r, double shp[]) {
    shp[0] = (1 - r) / 2;
    shp[1] = (1 + r) / 2;
  }

  void dShapeLine2(double r, double ds[]) {
    ds[0] = -.5;
    ds[1] = .5;
  }
}


PDE1dImpl::PDE1dImpl(PDE1dDefn &pde, PDE1dOptions &options) : 
pde(pde), options(options)
{
  ida = 0;
  intRule = new GausLegendreIntRule(3);
  mesh = pde.getMesh();
  checkIncreasing(mesh, 5, "meshPts");
  tspan = pde.getTimeSpan();
  checkIncreasing(tspan, 6, "timePts");
  numNodes = mesh.size();
  numTimes = tspan.size();
  numDepVars = pde.getNumEquations();
  numFEMEqns = numNodes*numDepVars;
  dirConsFlagsLeft.resize(numDepVars);
  dirConsFlagsRight.resize(numDepVars);
  y0.resize(numDepVars, numNodes);
  // get initial conditions
  RealVector ic(numDepVars);
  for (int i = 0; i < numNodes; i++) {
    pde.evalIC(mesh(i), ic);
    y0.col(i) = ic;
  }
  // flag dirichlet constraints
  int nm1 = numNodes - 1;
  RealVector ul = y0.col(0), ur = y0.col(nm1);
  bc.pl.resize(numDepVars);
  bc.pr.resize(numDepVars);
  bc.ql.resize(numDepVars);
  bc.qr.resize(numDepVars);
  pde.evalBC(mesh(0), ul, mesh(nm1), ur, tspan(0), bc);
  for (int i = 0; i < numDepVars; i++) {
    dirConsFlagsLeft[i] = bc.ql[i] == 0;
    dirConsFlagsRight[i] = bc.qr[i] == 0;
  }

  coeffs.c.resize(numDepVars);
  coeffs.f.resize(numDepVars);
  coeffs.s.resize(numDepVars);

}


PDE1dImpl::~PDE1dImpl()
{
  if(ida) IDAFree(&ida);
  delete intRule;
}

namespace {

  void print(N_Vector v, const char *title) {
    double *d = NV_DATA_S(v);
    int len = NV_LENGTH_S(v);
    int count = 0;
    printf("Vector: %s(%d)\n", title, len);
    for (int i = 0; i < len; i++) {
      printf(" %14.9e", d[i]);
      if (++count == 6) {
        printf("\n");
        count = 0;
      }
    }
    if (count)
      printf("\n");
  }

  int resFunc(realtype tres, N_Vector uu, N_Vector up, N_Vector resval,
    void *user_data) {

    PDE1dImpl *pde = (PDE1dImpl*) user_data;

    MapVec u(NV_DATA_S(uu), NV_LENGTH_S(uu));
    MapVec uDot(NV_DATA_S(up), NV_LENGTH_S(up));
    MapVec res(NV_DATA_S(resval), NV_LENGTH_S(resval));

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

    fem->calcJacobianODE(t, c_j, uu, up, r, Jac);
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

  template<class T>
  bool initConditionsChanged(const T &ic, const T &icMod, double tol) {
 
    for (int i = 0; i < ic.size(); i++) {
      if (std::abs(ic[i] - icMod[i]) > tol) {
        //printf("IC: initial=%g, modified=%g\n", ic[i], icMod[i]);
        return true;
      }
    }
    return false;
  }

}

int PDE1dImpl::solveTransient(PDESolution &sol)
{
  sol.u.resize(numTimes, numFEMEqns);
  sol.time.resize(numTimes);

  int neqImpl = numFEMEqns;
  /* Create vectors uu, up, res, constraints, id. */
  SunVector uu(neqImpl), up(neqImpl), res(neqImpl), id(neqImpl);
  double *udata = &uu[0];
  double *updata = &up[0];
  MapVec u(udata, neqImpl);
  std::copy_n(y0.data(), numFEMEqns, udata);
  up.setConstant(0);
#if 0
  // test for consistent initial conditions
  resFunc(0, uu(), up(), res(), this);
  MapVec resVec(&res[0], neqImpl);
  double initResNorm = resVec.dot(resVec);
  printf("initResNorm=%12.3e\n", sqrt(initResNorm));
#endif
  setAlgVarFlags(id());

  /* Call IDACreate and IDAMalloc to initialize solution */
  ida = IDACreate();
  check_flag(ida, "IDACreate", 0);

  int ier = IDASetUserData(ida, this);
  check_flag(&ier, "IDASetUserData", 1);
#if ! DAE_Y_INIT
  ier = IDASetId(ida, id());
  check_flag(&ier, "IDASetId", 1);
#endif
  double t0 = 0, tf = .05;
  IDAResFn resFn = resFunc;
  ier = IDAInit(ida, resFn, t0, uu(), up());
  check_flag(&ier, "IDAInit", 1);
  const double relTol = options.getRelTol(), absTol = options.getAbsTol();
  ier = IDASStolerances(ida, relTol, absTol);
  check_flag(&ier, "IDASStolerances", 1);
  ier = IDASetMaxNumSteps(ida, options.getMaxSteps());
  check_flag(&ier, "IDASetMaxNumSteps", 1);
#if SUN_USING_SPARSE
  //printf("Using sparse solver.\n");
  SparseMat P;
  calcJacPattern(P);
  //cout << P.toDense() << endl;
  ier = IDAKLU(ida, numFEMEqns, P.nonZeros());
  check_flag(&ier, "IDAKLU", 1);
  ier = IDASlsSetSparseJacFn(ida, jacFunc);
  check_flag(&ier, "IDASlsSetSparseJacFn", 1);
  fDiffJac = new FiniteDiffJacobian(P);
#elif BAND_SOLVER
  /* Call IDABand to specify the linear solver. */
  int mu = 2*numDepVars-1, ml = mu;
  ier = IDABand(ida, neqImpl, mu, ml);
  check_flag(&ier, "IDABand", 1);
#else
  ier = IDADense(ida, neqImpl);
  check_flag(&ier, "IDADense", 1);
#endif

  // Call IDACalcIC to correct the initial values. 
#if DAE_Y_INIT
  ier = IDACalcIC(ida, IDA_Y_INIT, tf);
#else
  ier = IDACalcIC(ida, IDA_YA_YDP_INIT, tf);
#endif
  if (ier < 0) {
    throw PDE1dException("pde1d:consistent_ic",
      "Unable to calculate a consistent initial solution to the PDE.\n"
      "Often this is caused by an incorrect specification of the boundary conditions.");
  }
  SunVector yy0_mod(neqImpl), yp0_mod(neqImpl);
  ier = IDAGetConsistentIC(ida, yy0_mod(), yp0_mod());
  check_flag(&ier, "IDAGetConsistentIC", 1);
#if 0
  SunVector diff_vec(neqImpl), unit_vec(neqImpl);
  N_VLinearSum(1, uu(), -1, yy0_mod(), diff_vec());
  unit_vec.setConstant(1);
  double normDiff = N_VWrmsNorm(diff_vec(), unit_vec());
  printf("unit_vec=%12.3e\n", normDiff);
#endif
  MapVec uMod(&yy0_mod[0], neqImpl);
  if (initConditionsChanged(u, uMod, 100*absTol)) {
    PDE1dWarningMsg("pde1d:init_cond_changed",
   "User-defined initial conditions were changed "
      "to create a consistent solution to the equations at "
      "the initial time.\n");
  }
#if 0
  print(uu(), "yy0");
  print(yy0_mod(), "yy0_mod");
  print(up(), "yp0");
  print(yp0_mod(), "yp0_mod");
  // check up
  //resFn(0, yy0_mod, yp0_mod, res, this);
  //print(res, "res");
#endif
  
  sol.time(0) = tspan(0);
  MapVec y0Mod(&yy0_mod[0], neqImpl);
  sol.u.row(0) = y0Mod;
  for (int i = 1; i < numTimes; i++) {
    double tout=tspan(i), tret;
    ier = IDASolve(ida, tout, &tret, uu(), up(), IDA_NORMAL);
    if (ier < 0) {
      printStats();
      char msg[1024];
      sprintf(msg, "Time integration failed at t=%15.6e before reaching final time.\n"
        "Often this is caused by one or more dependent variables becoming unbounded.",
        tret);
      throw PDE1dException("pde1d:integ_failure", msg);
    }
    sol.time(i) = tret;
    sol.u.row(i) = u;
  }
  printStats();

  return 0;
}

template<class T, class TR>
void PDE1dImpl::calcGlobalEqns(double t, T &u, T &up, 
  TR &Cxd, TR &F, TR &S)
{
  const int nen = numElemNodes; // work around for g++ link error
  const int m = pde.getCoordSystem();
  int numEqns = numFEMEqns;
  int numElemEqns = numDepVars*nen;
  Cxd.setZero();

  //cout << "u=" << u.transpose() << endl;
  MapMat u2(u.data(), numDepVars, numNodes);
  MapMat up2(up.data(), numDepVars, numNodes);
  //cout << "u2=" << u2.transpose() << endl;
  //cout << "up2=" << up2.transpose() << endl;

  int numIntPts = intRule->getNumPoints();
  RealVector intPts(numIntPts), intWts(numIntPts);
  RealMatrix N(nen, numIntPts), dN(nen, numIntPts);
  for (int i = 0; i < numIntPts; i++) {
    double pt;
    intRule->getPoint(i, pt, intWts[i]);
    intPts[i] = pt;
    shapeLine2(pt, N.col(i).data());
    dShapeLine2(pt, dN.col(i).data());
  }
#if 0
  cout << "intPts=" << intPts.transpose() << endl;
  cout << "intWts=" << intWts.transpose() << endl;
  cout << "N\n" << N << endl;
#endif

  RealMatrix eCMat = RealMatrix::Zero(numElemEqns, numElemEqns);
  RealVector eF = RealVector::Zero(numElemEqns);
  RealVector eS = RealVector::Zero(numElemEqns);
  RealVector eUp(numElemEqns);
  RealMatrix eCj(nen, nen);
  RealVector eFj(numElemNodes), eSj(numElemNodes);

  RealVector &x = mesh;
  int ne = numNodes - 1;
  int eInd = 0;
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
      pde.evalPDE(xi, t, ui, dUiDx, coeffs);
      checkCoeffs(coeffs);
      double xm = 1;
      if (m == 1)
        xm = xi;
      else if(m == 2)
        xm = xi*xi;
      coeffs.c *= xm;
      coeffs.s *= xm;
      coeffs.f *= xm;
 
      for (int j = 0; j < numDepVars; j++) {
        eCj  = N.col(i)*coeffs.c(j)*N.col(i).transpose()*jacWt;
        //eCMat += eCj;
        eFj = dNdx*coeffs.f(j)*jacWt;
        //eF += eFj;
        eSj = N.col(i)*coeffs.s(j)*jacWt;
        //eS += eSj;
        for (int k = 0; k < numElemNodes; k++) {
          int kk = j + k*numDepVars;
          eF(kk) += eFj(k);
          eS(kk) += eSj(k);
          for (int l = 0; l < numElemNodes; l++) {
            int ll = j + l*numDepVars;
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
    int numEqns = numFEMEqns;
    int numElemEqns = numDepVars*nen;
    Cxd.setZero();

    //cout << "u=" << u.transpose() << endl;
    MapMat u2(u.data(), numDepVars, numNodes);
    MapMat up2(up.data(), numDepVars, numNodes);
    //cout << "u2=" << u2.transpose() << endl;
    //cout << "up2=" << up2.transpose() << endl;

    int numIntPts = intRule->getNumPoints();
    RealVector intPts(numIntPts), intWts(numIntPts);
    RealMatrix N(nen, numIntPts), dN(nen, numIntPts);
    for (int i = 0; i < numIntPts; i++) {
      double pt;
      intRule->getPoint(i, pt, intWts[i]);
      intPts[i] = pt;
      shapeLine2(pt, N.col(i).data());
      dShapeLine2(pt, dN.col(i).data());
    }
#if 0
    cout << "intPts=" << intPts.transpose() << endl;
    cout << "intWts=" << intWts.transpose() << endl;
    cout << "N\n" << N << endl;
#endif

    // get the coefficients at all integ pts in a single call
    int ne = numNodes - 1;
    int numXPts = numIntPts*ne;
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
    coeffsAllPts.c.resize(numDepVars, numXPts);
    coeffsAllPts.f.resize(numDepVars, numXPts);
    coeffsAllPts.s.resize(numDepVars, numXPts);
    pde.evalPDE(xPts, t, uPts, duPts, coeffsAllPts);

    RealMatrix eCMat = RealMatrix::Zero(numElemEqns, numElemEqns);
    RealVector eF = RealVector::Zero(numElemEqns);
    RealVector eS = RealVector::Zero(numElemEqns);
    RealVector eUp(numElemEqns);
    RealMatrix eCj(nen, nen);
    RealVector eFj(numElemNodes), eSj(numElemNodes);

    
    ip = 0; 
    int eInd = 0;
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
        Eigen::Ref<Eigen::VectorXd> cIp = coeffsAllPts.c.col(ip);
        Eigen::Ref<Eigen::VectorXd> sIp = coeffsAllPts.s.col(ip);
        Eigen::Ref<Eigen::VectorXd> fIp = coeffsAllPts.f.col(ip);
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
            int kk = j + k*numDepVars;
            eF(kk) += eFj(k);
            eS(kk) += eSj(k);
            for (int l = 0; l < numElemNodes; l++) {
              int ll = j + l*numDepVars;
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


template<class T>
void PDE1dImpl::calcRHSODE(double time, T &u, T &up, T &R)
{
  RealVector Cxd = RealVector::Zero(numFEMEqns);
  RealVector F = RealVector::Zero(numFEMEqns);
  RealVector S = RealVector::Zero(numFEMEqns);

  if (options.isVectorized() && pde.hasVectorPDEEval())
    calcGlobalEqnsVec(time, u, up, Cxd, F, S);
  else
    calcGlobalEqns(time, u, up, Cxd, F, S);
  R = F - S;

  // apply constraints
  MapMat u2(u.data(), numDepVars, numNodes);
  int rightDofOff = numFEMEqns - numDepVars;
  int nm1 = numNodes - 1;
  RealVector ul = u2.col(0), ur = u2.col(nm1);
  pde.evalBC(mesh(0), ul, mesh(nm1), ur, time, bc);
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
  R += Cxd;
#if 0
  cout << "Cxd\n" << Cxd << endl;
  cout << "F\n" << F << endl;
#endif
}

void PDE1dImpl::testMats()
{
  RealVector u(numFEMEqns), up(numFEMEqns);
  for (int i = 0; i < numFEMEqns; i++) {
    u(i) = .1*(i+1);
    up(i) = 0;
  }
  RealVector Cxd = RealVector::Zero(numFEMEqns);
  RealVector F = RealVector::Zero(numFEMEqns);
  RealVector S = RealVector::Zero(numFEMEqns);
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


void PDE1dImpl::setAlgVarFlags(N_Vector id)
{
  N_VConst(1, id);
  double *iddata = NV_DATA_S(id);
  MapMat idMat(iddata, numDepVars, numNodes);

  // flag equations where c==0
  const int nen = numElemNodes;
  RealMatrix dN(nen, nen);
  double xi[] = { -1, 1 };
  for (int i = 0; i < numElemNodes; i++) {
    double pt = xi[i];
    dShapeLine2(pt, dN.col(i).data());
  }
 
  duPts.resize(numDepVars, numNodes);
  int nnm1 = numNodes - 1;
  for(int i = 0; i < numNodes; i++) {
    double L;
    if (i < nnm1)
      L = mesh(i + 1) - mesh(i);
    else
      L = mesh(i) - mesh(i - 1);
    double jac = L / 2;
    auto dNdx = dN.col(0) / jac;
    if (i < nnm1)
      duPts.col(i) = (y0.col(i)*dNdx(0) + y0.col(i + 1)*dNdx(1));
    else
      duPts.col(i) = (y0.col(i-1)*dNdx(0) + y0.col(i)*dNdx(1));
  }

  if (options.isVectorized() && pde.hasVectorPDEEval()) {
    coeffsAllPts.c.resize(numDepVars, numNodes);
    coeffsAllPts.f.resize(numDepVars, numNodes);
    coeffsAllPts.s.resize(numDepVars, numNodes);
    pde.evalPDE(mesh, 0, y0, duPts, coeffsAllPts);
    for (int i = 0; i < numNodes; i++) {
      for (int j = 0; j < numDepVars; j++)
        if (coeffsAllPts.c(j,i)==0)
          idMat(j, i) = 0;
    }
  }
  else {
    for (int i = 0; i < numNodes; i++) {
      double xi = mesh(i);
      const auto &ui = y0.col(i);
      const auto &dUiDx = duPts.col(i);
      pde.evalPDE(xi, 0, ui, dUiDx, coeffs);
      for (int j = 0; j < numDepVars; j++)
        if (coeffs.c[j] == 0)
          idMat(j, i) = 0;
    }
  }

  // account for dirichlet constraints at ends
  int rtBcOff = numFEMEqns - numDepVars;
  for (int i = 0; i < numDepVars; i++) {
    if (dirConsFlagsLeft[i])
      iddata[i] = 0;
    if (dirConsFlagsRight[i])
      iddata[i + rtBcOff] = 0;
  }
  //print(id, "id");
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

void PDE1dImpl::checkCoeffs(const PDE1dDefn::PDE &coeffs)
{
  for (int i = 0; i < coeffs.c.size(); i++) {
    if (coeffs.c[i] != 0)
      return;
  }
  throw PDE1dException("pde1d:no_parabolic_eqn", "At least one of the "
    "entries in the c-coefficient vector must be non-zero.");
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
  printf("Number of internal time steps = %d\n", nsteps);
  printf("Number of residual function calls = %d\n", nrevals);
  printf("Number of Jacobian calculations = %d\n", nlinsetups);
  printf("Number of solution accuracy test failures = %d\n", netfails);
  printf("Last internal time step size = %12.3e\n", hlast);
  printf("Number of nonlinear iterations = %d\n", nniters);
  printf("Number of nonlinear convergence failures = %d\n", nncfails);
}

void PDE1dImpl::calcJacPattern(Eigen::SparseMatrix<double> &J)
{
  J.resize(numFEMEqns, numFEMEqns);
  int n2 = numDepVars*numDepVars;
  int nel = numNodes - 1;
  int nnz = 3 * n2*(nel + 1); // approximate nnz
  J.reserve(nnz);
  int eOff = 0;
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
  J.makeCompressed();
}

template<class T, class T2>
void PDE1dImpl::calcJacobianODE(double time, double beta, T &u, T &up, T &res,
  T2 Jac) {
  fDiffJac->calcJacobian(time, 1, beta, u, up, res, resFunc, this, Jac);
}