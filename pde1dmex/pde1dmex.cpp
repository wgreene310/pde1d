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

#include <stdio.h>

#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

#include <mex.h>

#include "PDE1dMexInt.h"
#include "PDE1dImpl.h"
#include "PDE1dOptions.h"

typedef double Float;

namespace {
  template<class T>
  mxArray *toMxArray(const T &a)
  {
    mxArray *ma = mxCreateDoubleMatrix(a.rows(), a.cols(), mxREAL);
    double *dest = mxGetPr(ma);
    const double *src = a.data();
    std::copy_n(src, a.cols()*a.rows(), dest);
    return ma;
  }

  PDE1dOptions getOptions(const mxArray *opts) {
    if (!mxIsStruct(opts))
      mexErrMsgTxt("The 7th, options argument to " FUNC_NAME
      " must be a struct.");
    PDE1dOptions pdeOpts;
    int n = mxGetNumberOfFields(opts);
    for (int i = 0; i < n; i++) {
      const char *ni = mxGetFieldNameByNumber(opts, i);
      mxArray *val = mxGetFieldByNumber(opts, 0, i);
      if (boost::iequals(ni, "reltol")) {
        pdeOpts.setRelTol(mxGetScalar(val));
      }
      else if (boost::iequals(ni, "abstol")) {
        pdeOpts.setAbsTol(mxGetScalar(val));
      }
      else {
        char msg[1024];
        sprintf(msg, "The options argument contains the field \"%s\".\n"
          "This is not a currently-supported option and will be ignored.",
          ni);
        mexWarnMsgIdAndTxt("pdep1:unknown_option", msg);
      }
    }

    return pdeOpts;
  }
}

void mexFunction(int nlhs, mxArray*
  plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 6 && nrhs != 7)
    mexErrMsgTxt(FUNC_NAME " requires six or seven input arguments.");

  if (mxGetNumberOfElements(prhs[0]) != 1) {
    mexErrMsgTxt("First argument must be an integer scalar.");
  }

  int m = (int)mxGetScalar(prhs[0]);
  if (m != 0 && m != 1 && m != 2)
    mexErrMsgTxt("First argument must be either 0, 1, or 2");

  for (int i = 1; i < 4; i++) {
    if (!mxIsFunctionHandle(prhs[i])) {
      char msg[80];
      sprintf(msg, "Argument %d is not a function handle.", i + 1);
      mexErrMsgTxt(msg);
    }
  }

  if (mxGetNumberOfElements(prhs[4]) < 3)
    mexErrMsgTxt("Length of argument five, xmesh, must be at least three.");
  if (mxGetNumberOfElements(prhs[4]) < 3)
    mexErrMsgTxt("Length of argument siz, tspan, must be at least three.");

  PDE1dOptions opts;
  if (nrhs == 7)
    opts = getOptions(prhs[6]);

  if(nlhs != 1)
    mexErrMsgTxt("pde1d returns only a single matrix.");

  plhs[0] = 0;

  PDESolution pdeSol;
  int numPde = 0, numPts = 0;
  try {
    PDE1dMexInt pde(m, prhs[1], prhs[2], prhs[3],
      prhs[4], prhs[5]);
    numPde = pde.getNumEquations();
    numPts = pde.numNodes();
    PDE1dImpl pdeImpl(pde, opts);
    int err = pdeImpl.solveTransient(pdeSol);
    if (err)
      return;
  }
  catch (const std::exception &ex) {
    mexPrintf("Error in pde1d: %s.\n", ex.what());
    return;
  }
  catch (...) {
    mexPrintf("Internal error in pde1d.\n");
    return;
  }
  //mexPrintf("%d %d %d\n", pdeSol.time.size(), 
  //  pdeSol.u.rows(), pdeSol.u.cols());


  int numTimes = pdeSol.u.rows();
 /* mexPrintf("numTimes=%d, numPde=%d, numPts=%d\n",
    numTimes, numPde, numPts);*/
  mxArray *sol;
  if (numPde > 1) {
    int ndims, dims[] = { numTimes, numPts, numPde };
    ndims = 3;
    typedef boost::multi_array_ref<double, 3> Matrix3;
    boost::array<Matrix3::index, 3> shape =
    { { numTimes, numPde, numPts } };
    Matrix3 sol3(pdeSol.u.data(), shape, boost::fortran_storage_order());
    sol = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    double *p = mxGetPr(sol);
    for (int k = 0; k < numPde; k++) {
      for (int j = 0; j < numPts; j++) {
        for (int i = 0; i < numTimes; i++) {
          *p++ = sol3[i][k][j];
        }
      }
    }
  }
  else {
    sol = mxCreateDoubleMatrix(numTimes, numPts, mxREAL);
    std::copy_n(pdeSol.u.data(), numTimes*numPts*numPde, mxGetPr(sol));
  }

  plhs[0] = sol;

}