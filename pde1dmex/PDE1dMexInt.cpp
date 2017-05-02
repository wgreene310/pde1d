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

#include <iostream>
#include <stdio.h>

using std::cout;
using std::endl;

#include <boost/timer.hpp>

#include "PDE1dMexInt.h"
#include "MexInterface.h"

namespace {

  void print(const mxArray *a, const char *name)
  {
    size_t m = mxGetM(a);
    size_t n = mxGetN(a);
    const Eigen::Map<Eigen::MatrixXd> A(mxGetPr(a), m, n);
    mexPrintf("%s(%d,%d)\n", name, m, n);
    for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++)
        printf("%g ", A(i, j));
      printf("\n");
    }
  }

}

PDE1dMexInt::PDE1dMexInt(int m, const mxArray *pdefun, const mxArray *icfun,
  const mxArray *bcfun,
  const mxArray *xmesh, const mxArray *tspan) : mCoord(m), pdefun(pdefun),
  icfun(icfun), bcfun(bcfun), xmesh(xmesh), tspan(tspan)
{
  mxX1 = mxX2 = mxT = mxVec1 = mxVec2  = 0;
  size_t numMesh = mxGetNumberOfElements(xmesh);
  mesh.resize(numMesh);
  std::copy_n(mxGetPr(xmesh), numMesh, mesh.data());
  size_t numTime = mxGetNumberOfElements(tspan);

  tSpan.resize(numTime);
  std::copy_n(mxGetPr(tspan), numTime, tSpan.data());
  mxX1 = mxCreateDoubleScalar(0);
  mxX2 = mxCreateDoubleScalar(0);
  mxT = mxCreateDoubleScalar(0);
  setNumPde();
  mxVec1 = mxCreateDoubleMatrix(numPDE, 1, mxREAL);
  mxVec2 = mxCreateDoubleMatrix(numPDE, 1, mxREAL);
  odefun = 0;
  odeIcFun = 0;
  odeMesh = 0;
  numODE = 0;
  mxV = 0;
  mxVDot = 0;
  mxOdeU = 0;
  mxOdeDuDx = 0;
  mxOdeR = 0;
  mxOdeDuDt = 0;
  mxOdeDuDxDt = 0;
}


PDE1dMexInt::~PDE1dMexInt()
{
  mxDestroyArray(mxX1);
  mxDestroyArray(mxX2);
  mxDestroyArray(mxT);
  mxDestroyArray(mxVec1);
  mxDestroyArray(mxVec2);
  destroy(mxV);
  destroy(mxVDot);
  destroy(mxOdeU);
  destroy(mxOdeDuDx);
  destroy(mxOdeR);
  destroy(mxOdeDuDt);
  destroy(mxOdeDuDxDt);
}

void PDE1dMexInt::destroy(mxArray *a)
{
  // not sure if mex does this test?
  if (a)
    mxDestroyArray(a);
}

void PDE1dMexInt::setODEDefn(const mxArray *odeFun, const mxArray *icFun,
  const mxArray *odemesh)
{
  odefun = odeFun;
  odeIcFun = icFun;
  odeMesh = odemesh;
  setNumOde();
  size_t numXi = mxGetNumberOfElements(odeMesh);
  mxV = mxCreateDoubleMatrix(numODE, 1, mxREAL);
  mxVDot = mxCreateDoubleMatrix(numODE, 1, mxREAL);
  mxOdeU = mxCreateDoubleMatrix(numXi, 1, mxREAL);
  mxOdeDuDx = mxCreateDoubleMatrix(numXi, 1, mxREAL);
  mxOdeR = mxCreateDoubleMatrix(numXi, 1, mxREAL);
  mxOdeDuDt = mxCreateDoubleMatrix(numXi, 1, mxREAL);
  mxOdeDuDxDt = mxCreateDoubleMatrix(numXi, 1, mxREAL);
  odeMeshVec = MexInterface::fromMxArrayVec(odeMesh);
}

void PDE1dMexInt::evalIC(double x, RealVector &ic)
{
  setScalar(x, mxX1);
  const int nargout = 1, nargin = 2;
  const mxArray *funcInp[] = { icfun, mxX1 };
  RealVector *outArgs[] = { &ic };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalODEIC(RealVector &ic)
{
  const int nargout = 1, nargin = 1;
  const mxArray *funcInp[] = { odeIcFun };
  RealVector *outArgs[] = { &ic };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalBC(double xl, const RealVector &ul,
  double xr, const RealVector &ur, double t, 
  const RealVector &v, const RealVector &vDot, BC &bc) {
  setScalar(xl, mxX1);
  setVector(ul, mxVec1);
  setScalar(xr, mxX2);
  setVector(ur, mxVec2);
  setScalar(t, mxT);
  const int nargout = 4;
  int nargin = 6;
  if (numODE) {
    nargin = 8;
    setVector(v, mxV);
    setVector(vDot, mxVDot);
  }
  const mxArray *funcInp[] = { bcfun, mxX1, mxVec1, mxX2, mxVec2, mxT,
    mxV, mxVDot };
  RealVector *outArgs[] = { &bc.pl, &bc.ql, &bc.pr, &bc.qr };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalPDE(double x, double t,
  const RealVector &u, const RealVector &DuDx, 
  const RealVector &v, const RealVector &vDot, PDE &pde) {
  // Evaluate pde coefficients one point at a time
  // [c,f,s] = heatpde(x,t,u,DuDx)
  setScalar(x, mxX1);
  setScalar(t, mxT);
  setVector(u, mxVec1);
  setVector(DuDx, mxVec2);
  const int nargout = 3;
  int nargin = 5;
  if (numODE) {
    nargin = 7;
    setVector(v, mxV);
    setVector(vDot, mxVDot);
  }
  const mxArray *funcInp[] = { pdefun, mxX1, mxT, mxVec1, mxVec2,
    mxV, mxVDot };
  RealVector *outArgs[] = { &pde.c, &pde.f, &pde.s };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalPDE(const RealVector &x, double t,
  const RealMatrix &u, const RealMatrix &DuDx, 
  const RealVector &v, const RealVector &vDot, PDECoeff &pde)
{
  // Evaluate pde coefficients at all x-locations
  // [c,f,s] = heatpde(x,t,u,DuDx)
  setMatrix(x.transpose(), mxX1);
  setScalar(t, mxT);
  setMatrix(u, mxVec1);
  setMatrix(DuDx, mxVec2);
  const int nargout = 3;
  int nargin = 5;
  if (numODE) {
    nargin = 7;
    setVector(v, mxV);
    setVector(vDot, mxVDot);
  }
  const mxArray *funcInp[] = { pdefun, mxX1, mxT, mxVec1, mxVec2,
  mxV, mxVDot};
  RealMatrix *outArgs[] = { &pde.c, &pde.f, &pde.s };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalODE(double t, const RealVector &v,
  const RealVector &vdot, const RealMatrix &u, const RealMatrix &DuDx,
  const RealMatrix &odeR, const RealMatrix &odeDuDt,
  const RealMatrix &odeDuDxDt, RealVector &f)
{
  // odeFunc(t,v,vdot,x,u,DuDx)
  setScalar(t, mxT);
  setVector(v, mxV);
  setVector(vdot, mxVDot);
  setMatrix(u, mxOdeU);
  setMatrix(DuDx, mxOdeDuDx);
  setMatrix(odeR, mxOdeR);
  setMatrix(odeDuDt, mxOdeDuDt);
  setMatrix(odeDuDxDt, mxOdeDuDxDt);
  const int nargout = 1;
  const int nargin = 10;
  const mxArray *funcInp[] = { odefun, mxT, mxV, mxVDot, odeMesh,
    mxOdeU, mxOdeDuDx, mxOdeR, mxOdeDuDt, mxOdeDuDxDt };
  RealVector *outArgs[] = { &f };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

const RealVector &PDE1dMexInt::getMesh() const {
  return mesh;
}

const RealVector &PDE1dMexInt::getODEMesh()
{
  return odeMeshVec;
}

const RealVector &PDE1dMexInt::getTimeSpan() const { return tSpan; }

void PDE1dMexInt::callMatlab(const mxArray *inArgs[], int nargin,
  RealVector *outArgs[], int nargout)
{
  int err = mexCallMATLAB(nargout, matOutArgs, nargin,
    const_cast<mxArray**>(inArgs), "feval");
  if (err) {
    char msg[1024];
    std::string funcName = getFuncNameFromHandle(inArgs[0]);
    sprintf(msg, "An error occurred in the call to user-defined function:\n\"%s\".",
      funcName.c_str());
    mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:err", msg);
  }
  for (int i = 0; i < nargout; i++) {
    mxArray *a = matOutArgs[i];
    if (! a)
      mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:arg", "Error in mexCallMATLAB arg.");
    size_t retLen = mxGetNumberOfElements(a);
    size_t exLen = outArgs[i]->size();
    if (retLen != exLen) {
      char msg[1024];
      std::string funcName = getFuncNameFromHandle(inArgs[0]);
      size_t m = mxGetM(a), n = mxGetN(a);
      sprintf(msg, "In the call to user-defined function:\n\"%s\"\n"
        "returned entry %d had size (%zd x %zd) but a vector of size (%zd x 1)"
        " was expected.", funcName.c_str(), i + 1, m, n, exLen);
      mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:arglen", msg);
    }
    std::copy_n(mxGetPr(a), retLen, outArgs[i]->data());
    if (a)
      mxDestroyArray(a);
  }
}

void PDE1dMexInt::callMatlab(const mxArray *inArgs[], int nargin,
  RealMatrix *outArgs[], int nargout)
{
  int err = mexCallMATLAB(nargout, matOutArgs, nargin,
    const_cast<mxArray**>(inArgs), "feval");
  if (err) {
    char msg[1024];
    std::string funcName = getFuncNameFromHandle(inArgs[0]);
    sprintf(msg, "An error occurred in the call to user-defined function:\n\"%s\".",
      funcName.c_str());
    mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:err", msg);
  }
  for (int i = 0; i < nargout; i++) {
    mxArray *a = matOutArgs[i];
    if (!a)
      mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:arg", "Error in mexCallMATLAB arg.");
    size_t retRows = mxGetM(a);
    size_t retCols = mxGetN(a);
    size_t exRows = outArgs[i]->rows();
    size_t exCols = outArgs[i]->cols();
    if (retRows != exRows || retCols != exCols) {
      char msg[1024];
      std::string funcName = getFuncNameFromHandle(inArgs[0]);
      int m = mxGetM(a), n = mxGetN(a);
      sprintf(msg, "In the call to user-defined function:\n\"%s\"\n"
        "returned entry %d had size (%zd x %zd) but a matrix of size (%zd x %zd)"
        " was expected.", funcName.c_str(), i + 1, retRows, retCols, 
        exRows, exCols);
      mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:arglen", msg);
    }
    std::copy_n(mxGetPr(a), retRows*retCols, outArgs[i]->data());
    if (a)
      mxDestroyArray(a);
  }
}

void PDE1dMexInt::callMatlab(const mxArray *inArgs[], int nargin, int nargout)
{
  int err = mexCallMATLAB(nargout, matOutArgs, nargin,
    const_cast<mxArray**>(inArgs), "feval");
  if (err) {
    char msg[1024];
    std::string funcName = getFuncNameFromHandle(inArgs[0]);
    sprintf(msg, "An error occurred in the call to user-defined function:\n\"%s\".",
      funcName.c_str());
    mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:err", msg);
  }
  for (int i = 0; i < nargout; i++) {
    mxArray *a = matOutArgs[i];
    if (!a)
      mexErrMsgIdAndTxt("pde1d:mexCallMATLAB:arg", "Error in mexCallMATLAB arg.");
  }
}

void PDE1dMexInt::setNumPde() {
  setScalar(mesh(0), mxX1);
  mxArray *initCond = 0;
  const mxArray *funcInp[] = { icfun, mxX1 };
  int err = mexCallMATLAB(1, &initCond, 2,
    const_cast<mxArray**>(funcInp), "feval");
  if (err)
    mexErrMsgIdAndTxt("pde1d:mexCallMATLAB",
    "Error in mexCallMATLAB.\n");
  numPDE = mxGetNumberOfElements(initCond);
  if (initCond)
    mxDestroyArray(initCond);
}

void PDE1dMexInt::setNumOde()
{
  mxArray *initCond = 0;
  const mxArray *funcInp[] = { odeIcFun};
  int err = mexCallMATLAB(1, &initCond, 1,
    const_cast<mxArray**>(funcInp), "feval");
  if (err)
    mexErrMsgIdAndTxt("pde1d:mexCallMATLAB",
    "Error in mexCallMATLAB.\n");
  numODE = mxGetNumberOfElements(initCond);
  if (initCond)
    mxDestroyArray(initCond);
  //printf("numODE=%d\n", numODE);
}

std::string PDE1dMexInt::getFuncNameFromHandle(const mxArray *fh)
{
  mxArray *funcName = 0;
  int err = mexCallMATLAB(1, &funcName, 1,
    const_cast<mxArray**>(&fh), "func2str");
  if (err)
    mexErrMsgIdAndTxt("pde1d:mexCallMATLAB",
    "Error in mexCallMATLAB.\n");
  const int bufLen = 1024;
  char buf[bufLen];
  int len = mxGetString(funcName, buf, bufLen);
  std::string nam(buf);
  return nam;
}