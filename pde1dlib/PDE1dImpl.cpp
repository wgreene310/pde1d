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

#include <ida/ida.h>
#include <ida/ida_band.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

#include "PDE1dImpl.h"
#include "PDE1dOptions.h"

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

#define DEBUG_MATS 0

namespace {

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
  intRule = new GausLegendreIntRule(3);
  mesh = pde.getMesh();
  tspan = pde.getTimeSpan();
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
  delete intRule;
}

#define BAND_SOLVER 1

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
    printf("heatres called\n");
    print(uu, "u");
    print(up, "up");
#endif

    pde->calcRHSODE(tres, u, uDot, res);


#if DEBUG
    print(resval, "resval");
    exit(1);
#endif

    return 0;
  }

  int check_flag(void *flagvalue, const char *funcname, int opt)
  {
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
      fprintf(stderr,
        "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      return(1);
    }
    else if (opt == 1) {
      /* Check if flag < 0 */
      errflag = (int *)flagvalue;
      if (*errflag < 0) {
        fprintf(stderr,
          "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
        return(1);
      }
    }
    else if (opt == 2 && flagvalue == NULL) {
      /* Check if function returned NULL pointer - no memory allocated */
      fprintf(stderr,
        "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      return(1);
    }

    return(0);
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
  N_Vector uu = N_VNew_Serial(neqImpl);
  if (check_flag((void *)uu, "N_VNew_Serial", 0)) return(1);
  N_Vector up = N_VNew_Serial(neqImpl);
  if (check_flag((void *)up, "N_VNew_Serial", 0)) return(1);
  N_Vector res = N_VNew_Serial(neqImpl);
  if (check_flag((void *)res, "N_VNew_Serial", 0)) return(1);
  N_Vector id = N_VNew_Serial(neqImpl);
  if (check_flag((void *)id, "N_VNew_Serial", 0)) return(1);
  double *udata = NV_DATA_S(uu);
  double *updata = NV_DATA_S(up);
  MapVec u(udata, neqImpl);
  std::copy_n(y0.data(), numFEMEqns, udata);
  N_VConst(0, up);
  setAlgVarFlags(id);

  /* Call IDACreate and IDAMalloc to initialize solution */
  void *ida = IDACreate();
  if (check_flag((void *)ida, "IDACreate", 0)) return(1);

  int ier = IDASetUserData(ida, this);
  if (check_flag(&ier, "IDASetUserData", 1)) return(1);
  ier = IDASetId(ida, id);
  if (check_flag(&ier, "IDASetId", 1)) return(1);
  double t0 = 0, tf = .05;
  IDAResFn resFn = resFunc;
  ier = IDAInit(ida, resFn, t0, uu, up);
  if (check_flag(&ier, "IDAInit", 1)) return(1);
  const double relTol = options.getRelTol(), absTol = options.getAbsTol();
  ier = IDASStolerances(ida, relTol, absTol);
  if (check_flag(&ier, "IDASStolerances", 1)) return(1);
  /* Call IDABand to specify the linear solver. */
  int mu = numDepVars, ml = numDepVars;
  ier = IDABand(ida, neqImpl, mu, ml);
  if (check_flag(&ier, "IDABand", 1)) return(1);
#if 1
  /* Call IDACalcIC to correct the initial values. */
  ier = IDACalcIC(ida, IDA_YA_YDP_INIT, tf);
  if (check_flag(&ier, "IDACalcIC", 1)) return(1);
  N_Vector yy0_mod = N_VNew_Serial(neqImpl);
  N_Vector yp0_mod = N_VNew_Serial(neqImpl);
  ier = IDAGetConsistentIC(ida, yy0_mod, yp0_mod);
  if (check_flag(&ier, "IDAGetConsistentIC", 1)) return(1);
  MapVec uMod(NV_DATA_S(yy0_mod), neqImpl);
  if (initConditionsChanged(u, uMod, absTol)) {
    printf("User-defined initial conditions were changed "
      "to create a consistent solution to the equations at "
      "the initial time.\n");
  }
#if 0
  print(uu, "yy0");
  print(yy0_mod, "yy0_mod");
  // check up
  //resFn(0, yy0_mod, yp0_mod, res, this);
  //print(res, "res");
#endif

#endif
  
  sol.time(0) = tspan(0);
  MapVec y0Mod(NV_DATA_S(yy0_mod), neqImpl);
  sol.u.row(0) = y0Mod;
  for (int i = 1; i < numTimes; i++) {
    double tout=tspan(i), tret;
    ier = IDASolve(ida, tout, &tret, uu, up, IDA_NORMAL);
    if (check_flag(&ier, "IDASolve", 1)) return(1);
    sol.time(i) = tret;
    sol.u.row(i) = u;
  }

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
    for (int i = 0; i < numIntPts; i++) {
      // assume two node elements
      double jacWt = jac*intWts(i);
      double xi = x(e)*N(0, i) + x(e + 1)*N(1, i);
      auto ui = u2.col(e)*N(0, i) + u2.col(e + 1)*N(1, i);
      auto dNdx = dN.col(0) / jac;
      auto dUiDx = (u2.col(e)*dNdx(0) + u2.col(e + 1)*dNdx(1));
      pde.evalPDE(xi, t, ui, dUiDx, coeffs);
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
#if 0
        int off1 = j;
        for (int k = 0; k < numElemNodes; k++) {
          eF(off1 + k) += eFj(k);
          eS(off1 + k) += eSj(k);
          int off2 = j;
          for (int l = 0; l < numElemNodes; l++) {
            eCMat(off1 + k, off2 + l) += eCj(k, l);
            off2 += numDepVars-1;
          }
          off1 += numDepVars-1;
        }
#else
        for (int k = 0; k < numElemNodes; k++) {
          int kk = j + k*numDepVars;
          eF(kk) += eFj(k);
          eS(kk) += eSj(k);
          for (int l = 0; l < numElemNodes; l++) {
            int ll = j + l*numDepVars;
            eCMat(kk, ll) += eCj(k, l);
          }
        }
#endif
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

template<class T>
void PDE1dImpl::calcRHSODE(double time, T &u, T &up, T &R)
{
  RealVector Cxd = RealVector::Zero(numFEMEqns);
  RealVector F = RealVector::Zero(numFEMEqns);
  RealVector S = RealVector::Zero(numFEMEqns);

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
    u(i) = i;
    up(i) = i;
  }
  RealVector Cxd = RealVector::Zero(numFEMEqns);
  RealVector F = RealVector::Zero(numFEMEqns);
  RealVector S = RealVector::Zero(numFEMEqns);
  calcGlobalEqns(0, u, up, Cxd, F, S);
  cout << "Cxd\n" << Cxd.transpose() << endl;
  cout << "F\n" << F.transpose() << endl;
  cout << "S\n" << S.transpose() << endl;
  RealVector R = RealVector::Zero(numFEMEqns);
  calcRHSODE(0, u, up, R);
  cout << "R\n" << R.transpose() << endl;
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
  RealVector dUiDx(numDepVars);
  int nnm1 = numNodes - 1;
  for(int i = 0; i < numNodes; i++) {
    double xi = mesh(i);
    auto ui = y0.col(i);
    double L;
    if (i < nnm1)
      L = mesh(i + 1) - mesh(i);
    else
      L = mesh(i) - mesh(i - 1);
    double jac = L / 2;
    auto dNdx = dN.col(0) / jac;
    if (i < nnm1)
      dUiDx = (y0.col(i)*dNdx(0) + y0.col(i + 1)*dNdx(1));
    else
      dUiDx = (y0.col(i-1)*dNdx(0) + y0.col(i)*dNdx(1));
    pde.evalPDE(xi, 0, ui, dUiDx, coeffs);
    for (int j = 0; j < numDepVars; j++)
      if (coeffs.c[j] == 0)
        idMat(j, i) = 0;
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