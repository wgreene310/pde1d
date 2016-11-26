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
#include <stdexcept>

#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

#include <mex.h>

#include "MexInterface.h"
#include "PDE1dMexInt.h"
#include "PDE1dImpl.h"
#include "PDE1dOptions.h"
#include "PDE1dException.h"
#include "PDE1dWarningMsg.h"

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
      mexErrMsgIdAndTxt("pde1d:optins_type", 
      "The last options argument to " FUNC_NAME " must be a struct.");
    PDE1dOptions pdeOpts;
    int n = mxGetNumberOfFields(opts);
    for (int i = 0; i < n; i++) {
      const char *ni = mxGetFieldNameByNumber(opts, i);
      mxArray *val = mxGetFieldByNumber(opts, 0, i);
      if (boost::iequals(ni, "reltol")) 
        pdeOpts.setRelTol(mxGetScalar(val));
      else if (boost::iequals(ni, "abstol"))
        pdeOpts.setAbsTol(mxGetScalar(val));
      else if (boost::iequals(ni, "vectorized")) {
        const int buflen = 1024;
        char buf[buflen];
        mxGetString(val, buf, buflen);
        bool isVec;
        if (boost::iequals(buf, "on"))
          isVec = true;
        else if (boost::iequals(buf, "off"))
          isVec = false;
        else
          mexErrMsgIdAndTxt("pde1d:invalidVectorized",
          "The value of the \"Vectorized\" option must be either \"On\" or \"Off\".");
        pdeOpts.setVectorized(isVec);
      }
      else if (boost::iequals(ni, "maxsteps")) {
        int mxs = (int)mxGetScalar(val);
        pdeOpts.setMaxSteps(mxs);
      }
      else if (boost::iequals(ni, "stats")) {
        const int buflen = 1024;
        char buf[buflen];
        mxGetString(val, buf, buflen);
        bool doStats;
        if (boost::iequals(buf, "on"))
          doStats = true;
        else if (boost::iequals(buf, "off"))
          doStats = false;
        else
          mexErrMsgIdAndTxt("pde1d:invalidStats",
          "The value of the \"Stats\" option must be either \"On\" or \"Off\".");
        pdeOpts.setPrintStats(doStats);
      }
      else if (boost::iequals(ni, "icmethod")) {
        int icMethod = (int) mxGetScalar(val);
        pdeOpts.setICMethod(icMethod);
      }
      else if (boost::iequals(ni, "icdiagnostics")) {
        int icDiag = (int) mxGetScalar(val);
        pdeOpts.setICDiagnostics(icDiag);
      }
      else {
        char msg[1024];
        sprintf(msg, "The options argument contains the field \"%s\".\n"
          "This is not a currently-supported option and will be ignored.",
          ni);
        mexWarnMsgIdAndTxt("pde1d:unknown_option", msg);
      }
    }

    return pdeOpts;
  }
}

void PDE1dWarningMsg(const char *id, const char *msg) {
  mexWarnMsgIdAndTxt(id, msg);
}

/*
   solution = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts)
   solution = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts,options)
   solution = pde1d(m,pdeFunc,icFunc,bcFunc,meshPts,timePts,
      odefun,odeIcFunc,odeMesh)
*/

void mexFunction(int nlhs, mxArray*
  plhs[], int nrhs, const mxArray *prhs[])
{
  PDESolution pdeSol;
  int numPde = 0, numPts = 0, numOde = 0;
  try {

    //printf("nlhs=%d, nrhs=%d\n", nlhs, nrhs); return;
    int optsArg = -1;
    if (nrhs == 7)
      optsArg = 6;
    else if (nrhs == 10)
      optsArg = 9;
    else if (nrhs != 6 && nrhs != 9)
      mexErrMsgIdAndTxt("pde1d:nrhs",
      "Illegal number of input arguments passed to " FUNC_NAME);

    const bool hasODE = nrhs > 7;

    const mxArray *pM = prhs[0];
    if (!mxIsNumeric(pM) || mxGetNumberOfElements(pM) != 1) {
      mexErrMsgIdAndTxt("pde1d:invalid_m_type",
        "First argument must be an integer scalar.");
    }

    int m = (int)mxGetScalar(pM);
    if (m != 0 && m != 1 && m != 2)
      mexErrMsgIdAndTxt("pde1d:invalid_m_val",
      "First argument must be either 0, 1, or 2");

    for (int i = 1; i < 4; i++) {
      if (!mxIsFunctionHandle(prhs[i])) {
        char msg[80];
        sprintf(msg, "Argument %d is not a function handle.", i + 1);
        mexErrMsgIdAndTxt("pde1d:arg_not_func", msg);
      }
    }

    const mxArray *pX = prhs[4], *pT = prhs[5];
    if (!mxIsNumeric(pX) || mxIsComplex(pX))
      mexErrMsgIdAndTxt("pde1d:mesh_type",
      "Argument \"meshPts\" must be a real vector.");
    if (mxGetNumberOfElements(pX) < 3)
      mexErrMsgIdAndTxt("pde1d:mesh_length",
      "Length of argument \"meshPts\", must be at least three.");

    if (!mxIsNumeric(pT) || mxIsComplex(pT))
      mexErrMsgIdAndTxt("pde1d:time_type",
      "Argument \"timePts\" must be a real vector.");
    if (mxGetNumberOfElements(pT) < 3)
      mexErrMsgIdAndTxt("pde1d:time_length",
      "Length of argument \"timePts\", must be at least three.");

    PDE1dOptions opts;
    if (optsArg > 0)
      opts = getOptions(prhs[optsArg]);

    std::fill_n(plhs, nlhs, nullptr);

    PDE1dMexInt pde(m, prhs[1], prhs[2], prhs[3],
      prhs[4], prhs[5]);
    if (hasODE) {
      pde.setODEDefn(prhs[6], prhs[7], prhs[8]);
      if (nlhs > 2)
        mexErrMsgIdAndTxt("pde1d:nlhs",
        "pde1d returns only a two matrices when there are included ODEs.");
    }
    else {
      if (nlhs > 1)
        mexErrMsgIdAndTxt("pde1d:nlhs",
        "pde1d returns only a single matrix when there are no ODEs.");
    }
    numPde = pde.getNumPDE();
    numPts = pde.numNodes();
    numOde = pde.getNumODE();
    PDE1dImpl pdeImpl(pde, opts);
    int err = pdeImpl.solveTransient(pdeSol);
    if (err)
      return;
  }
  catch (const PDE1dException &ex) {
    mexErrMsgIdAndTxt(ex.getId(), ex.what());
  }
  catch (const std::exception &ex) {
    mexErrMsgIdAndTxt("pde1d:exception", ex.what());
  }
  catch (...) {
    mexErrMsgIdAndTxt("pde1d:internal_err", "Internal error in pde1d.\n");
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
  // odes are included
  if (nlhs == 2)
    plhs[1] = MexInterface::toMxArray(pdeSol.uOde);

}