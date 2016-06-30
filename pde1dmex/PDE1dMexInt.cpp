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

using std::cout;
using std::endl;

#include <boost/timer.hpp>

#include "PDE1dMexInt.h"

namespace {

  void print(const mxArray *a, const char *name)
  {
    int m = mxGetM(a);
    int n = mxGetN(a);
    const Eigen::Map<Eigen::MatrixXd> A(mxGetPr(a), m, n);
    mexPrintf("%s(%d,%d)\n", name, m, n);
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++)
        printf("%g ", A(i, j));
      printf("\n");
    }
  }

  template<class T>
  void print(const T &a, const char *name) {

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
}


PDE1dMexInt::~PDE1dMexInt()
{
  mxDestroyArray(mxX1);
  mxDestroyArray(mxX2);
  mxDestroyArray(mxT);
  mxDestroyArray(mxVec1);
  mxDestroyArray(mxVec2);
}

int PDE1dMexInt::getNumEquations()
{
  return numPDE;
}

void PDE1dMexInt::evalIC(double x, RealVector &ic)
{
  setScalar(x, mxX1);
  const int nargout = 1, nargin = 2;
  const mxArray *funcInp[] = { icfun, mxX1 };
  RealVector *outArgs[] = { &ic };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalBC(double xl, const RealVector &ul,
  double xr, const RealVector &ur, double t, BC &bc) {
  setScalar(xl, mxX1);
  setVector(ul, mxVec1);
  setScalar(xr, mxX2);
  setVector(ur, mxVec2);
  setScalar(t, mxT);
  const int nargout = 4, nargin = 6;
  const mxArray *funcInp[] = { bcfun, mxX1, mxVec1, mxX2, mxVec2, mxT };
  RealVector *outArgs[] = { &bc.pl, &bc.ql, &bc.pr, &bc.qr };
  callMatlab(funcInp, nargin, outArgs, nargout);
}

void PDE1dMexInt::evalPDE(double x, double t,
  const RealVector &u, const RealVector &DuDx, PDE &pde) {
  // [c,f,s] = heatpde(x,t,u,DuDx)
  setScalar(x, mxX1);
  setScalar(t, mxT);
  setVector(u, mxVec1);
  setVector(DuDx, mxVec2);
  const int nargout = 3, nargin = 5;
  const mxArray *funcInp[] = { pdefun, mxX1, mxT, mxVec1, mxVec2 };
  RealVector *outArgs[] = { &pde.c, &pde.f, &pde.s };
  callMatlab(funcInp, nargin, outArgs, nargout);
}
RealVector PDE1dMexInt::getMesh() {
  return mesh;
}

RealVector PDE1dMexInt::getTimeSpan() { return tSpan; }

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
    mexErrMsgIdAndTxt("MATLAB:callMatlab:err", msg);
  }
  for (int i = 0; i < nargout; i++) {
    mxArray *a = matOutArgs[i];
    if (! a)
      mexErrMsgIdAndTxt("MATLAB:callMatlab:arg", "Error in mexCallMATLAB arg.");
    int retLen = mxGetNumberOfElements(a);
    int exLen = outArgs[i]->size();
    if (retLen != exLen) {
      char msg[1024];
      std::string funcName = getFuncNameFromHandle(inArgs[0]);
      int m = mxGetM(a), n = mxGetN(a);
      sprintf(msg, "In the call to user-defined function:\n\"%s\"\n"
        "returned entry %d had size (%d x %d) but a vector of size (%d x 1)"
        " was expected.", funcName.c_str(), i + 1, m, n, exLen);
      mexErrMsgIdAndTxt("MATLAB:callMatlab:arglen", msg);
    }
    std::copy_n(mxGetPr(a), retLen, outArgs[i]->data());
    if (a)
      mxDestroyArray(a);
  }
}

void PDE1dMexInt::setNumPde() {
  setScalar(mesh(0), mxX1);
  mxArray *initCond = 0;
  const mxArray *funcInp[] = { icfun, mxX1 };
  int err = mexCallMATLAB(1, &initCond, 2,
    const_cast<mxArray**>(funcInp), "feval");
  if (err)
    mexErrMsgIdAndTxt("MATLAB:filterTriangles:minrhs",
    "Error in mexCallMATLAB.\n");
  numPDE = mxGetNumberOfElements(initCond);
  if (initCond)
    mxDestroyArray(initCond);
}

std::string PDE1dMexInt::getFuncNameFromHandle(const mxArray *fh)
{
  mxArray *funcName = 0;
  int err = mexCallMATLAB(1, &funcName, 1,
    const_cast<mxArray**>(&fh), "func2str");
  if (err)
    mexErrMsgIdAndTxt("MATLAB:filterTriangles:minrhs",
    "Error in mexCallMATLAB.\n");
  const int bufLen = 1024;
  char buf[bufLen];
  int len = mxGetString(funcName, buf, bufLen);
  std::string nam(buf);
  return nam;
}