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

#include "MexInterface.h"

MexInterface::MexInterface(int maxOutArgs) : maxOutArgs(maxOutArgs)
{
  matOutArgs.resize(maxOutArgs);
}


MexInterface::~MexInterface()
{
}

void MexInterface::print(const mxArray *a, const char *name)
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

void MexInterface::callMatlab(const mxArray *inArgs[], int nargin,
  RealVector *outArgs[], int nargout)
{
  int err = mexCallMATLAB(nargout, &matOutArgs[0], nargin,
    const_cast<mxArray**>(inArgs), "feval");
  if (err) {
    char msg[1024];
    std::string funcName = getFuncNameFromHandle(inArgs[0]);
    sprintf(msg, "An error occurred in the call to user-defined function:\n\"%s\".",
      funcName.c_str());
    mexErrMsgIdAndTxt("bvp1d:mexCallMATLAB:err", msg);
  }
  for (int i = 0; i < nargout; i++) {
    mxArray *a = matOutArgs[i];
    if (!a)
      mexErrMsgIdAndTxt("bvp1d:mexCallMATLAB:arg", "Error in mexCallMATLAB arg.");
    int retLen = mxGetNumberOfElements(a);
    int exLen = outArgs[i]->size();
    if (retLen != exLen) {
      char msg[1024];
      std::string funcName = getFuncNameFromHandle(inArgs[0]);
      int m = mxGetM(a), n = mxGetN(a);
      sprintf(msg, "In the call to user-defined function:\n\"%s\"\n"
        "returned entry %d had size (%d x %d) but a vector of size (%d x 1)"
        " was expected.", funcName.c_str(), i + 1, m, n, exLen);
      mexErrMsgIdAndTxt("bvp1d:mexCallMATLAB:arglen", msg);
    }
    std::copy_n(mxGetPr(a), retLen, outArgs[i]->data());
    if (a)
      mxDestroyArray(a);
  }
}

std::string MexInterface::getFuncNameFromHandle(const mxArray *fh)
{
  mxArray *funcName = 0;
  int err = mexCallMATLAB(1, &funcName, 1,
    const_cast<mxArray**>(&fh), "func2str");
  if (err)
    mexErrMsgIdAndTxt("bvp1d:mexCallMATLAB",
    "Error in mexCallMATLAB.\n");
  const int bufLen = 1024;
  char buf[bufLen];
  int len = mxGetString(funcName, buf, bufLen);
  std::string nam(buf);
  return nam;
}

Eigen::MatrixXd MexInterface::fromMxArray(const mxArray *a)
{
  const int m = mxGetM(a);
  const int n = mxGetN(a);
  Eigen::MatrixXd mat(m, n);
  std::copy_n(mxGetPr(a), m*n, mat.data());
  return mat;
}

Eigen::VectorXd MexInterface::fromMxArrayVec(const mxArray *a)
{
  const int m = mxGetM(a);
  const int n = mxGetN(a);
  Eigen::VectorXd mat(m*n);
  std::copy_n(mxGetPr(a), m*n, mat.data());
  return mat;
}