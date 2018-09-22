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

#include <stdio.h>
#include <stdexcept>
#include <iostream>
using std::cout;
using std::endl;

#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

#include <mex.h>

#include "MexInterface.h"
#include "PDE1dMexInt.h"
#include "PDE1dImpl.h"
#include "PDE1dOptions.h"
#include "PDESolution.h"
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

  void getOptions(const mxArray *opts, PDE1dOptions &pdeOpts,
    mxArray* &eventFunc) {
    if (!mxIsStruct(opts))
      pdeErrMsgIdAndTxt("pde1d:optins_type", 
      "The last options argument to " FUNC_NAME " must be a struct.");
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
          pdeErrMsgIdAndTxt("pde1d:invalidVectorized",
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
          pdeErrMsgIdAndTxt("pde1d:invalidStats",
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
      else if (boost::iequals(ni, "jacdiagnostics")) {
        int jacDiag = (int) mxGetScalar(val);
        pdeOpts.setJacDiagnostics(jacDiag);
      }
      else if (boost::iequals(ni, "polyorder")) {
        int porder = (int)mxGetScalar(val);
        pdeOpts.setPolyOrder(porder);
      }
      else if (boost::iequals(ni, "viewmesh")) {
        int vumesh = (int)mxGetScalar(val);
        pdeOpts.setViewMesh(vumesh);
      }
      else if (boost::iequals(ni, "diagonalMassMatrix")) {
        const int buflen = 1024;
        char buf[buflen];
        mxGetString(val, buf, buflen);
        bool useDiagMassMat;
        if (boost::iequals(buf, "on"))
          useDiagMassMat = true;
        else if (boost::iequals(buf, "off"))
          useDiagMassMat = false;
        else
          pdeErrMsgIdAndTxt("pde1d:invalidMassMat",
            "The value of the \"diagonalMassMatrix\" option must be either \"On\" or \"Off\".");
        pdeOpts.setDiagMassMat(useDiagMassMat);
      }
      else if (boost::iequals(ni, "events")) {
        if (!mxIsFunctionHandle(val))
          pdeErrMsgIdAndTxt("pde1d:invalidEventsFunc",
            "The value of the \"Events\" option must be a function handle.");
        eventFunc = val;
      }
      else {
        char msg[1024];
        sprintf(msg, "The options argument contains the field \"%s\".\n"
          "This is not a currently-supported option and will be ignored.",
          ni);
        mexWarnMsgIdAndTxt("pde1d:unknown_option", msg);
      }
    }
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
  try {
    //printf("nlhs=%d, nrhs=%d\n", nlhs, nrhs); return;
    // options struct is always the last argument in the
    // function call
    int optsArg = -1;
    if (nrhs == 7)
      optsArg = 6;
    else if (nrhs == 10)
      optsArg = 9;
    else if (nrhs != 6 && nrhs != 9)
      pdeErrMsgIdAndTxt("pde1d:nrhs",
        "Illegal number of input arguments passed to " FUNC_NAME);

    const bool hasODE = nrhs > 7;

    const mxArray *pM = prhs[0];
    if (!mxIsNumeric(pM) || mxGetNumberOfElements(pM) != 1) {
      pdeErrMsgIdAndTxt("pde1d:invalid_m_type",
        "First argument must be an integer scalar.");
    }

    int m = (int)mxGetScalar(pM);
    if (m != 0 && m != 1 && m != 2)
      pdeErrMsgIdAndTxt("pde1d:invalid_m_val",
      "First argument must be either 0, 1, or 2");

    for (int i = 1; i < 4; i++) {
      if (!mxIsFunctionHandle(prhs[i])) {
        char msg[80];
        sprintf(msg, "Argument %d is not a function handle.", i + 1);
        pdeErrMsgIdAndTxt("pde1d:arg_not_func", msg);
      }
    }

    const mxArray *pX = prhs[4];
    const mxArray *pT = prhs[5];
    if (!mxIsNumeric(pX) || mxIsComplex(pX))
      pdeErrMsgIdAndTxt("pde1d:mesh_type",
      "Argument \"meshPts\" must be a real vector.");
    if (mxGetNumberOfElements(pX) < 2)
      pdeErrMsgIdAndTxt("pde1d:mesh_length",
      "Length of argument \"meshPts\", must be at least two.");

    if (!mxIsNumeric(pT) || mxIsComplex(pT))
      pdeErrMsgIdAndTxt("pde1d:time_type",
      "Argument \"timePts\" must be a real vector.");
    if (mxGetNumberOfElements(pT) < 3)
      pdeErrMsgIdAndTxt("pde1d:time_length",
      "Length of argument \"timePts\", must be at least three.");

    PDE1dOptions opts;
    mxArray *eventsFunc = 0;
    if (optsArg > 0)
      getOptions(prhs[optsArg], opts, eventsFunc);

    std::fill_n(plhs, nlhs, nullptr);

    PDE1dMexInt pde(m, prhs[1], prhs[2], prhs[3],
      prhs[4], prhs[5]);
    pde.setEventsFunction(eventsFunc);
    if (hasODE) {
      pde.setODEDefn(prhs[6], prhs[7], prhs[8]);
      if (eventsFunc) {
        if (nlhs > 6)
          pdeErrMsgIdAndTxt("pde1d:nlhs",
            "pde1d returns six or fewer matrices when "
            "there are events and included ODEs.");
      }
      else {
        if (nlhs > 2)
          pdeErrMsgIdAndTxt("pde1d:nlhs",
            "pde1d returns only two matrices when there are included ODEs.");
      }
    }
    else {
      // no ODE
      if (eventsFunc) {
        if (nlhs > 5)
          pdeErrMsgIdAndTxt("pde1d:nlhs",
            "pde1d returns five or fewer matrices when "
            "there are events but no ODEs.");
      }
      else {
        if (nlhs > 1)
          pdeErrMsgIdAndTxt("pde1d:nlhs",
            "pde1d returns only a single matrix when there are no ODEs.");
      }
    }
    int numPde = pde.getNumPDE();
    int numOde = pde.getNumODE();
    PDE1dImpl pdeImpl(pde, opts);
    int viewMesh = opts.getViewMesh();
    PDESolution pdeSol(pde, pdeImpl.getModel(), viewMesh);
    int err = pdeImpl.solveTransient(pdeSol);
    if (err)
      return;

    int numTimes = pdeSol.numTimePoints();
    int numPts = pdeSol.numSpatialPoints();
#if 0
    printf("numTimes=%d, numPde=%d, numPts=%d\n",
    numTimes, numPde, numPts);
    pdeSol.print();
#endif
    mxArray *sol;
    if (numPde > 1) {
      mwSize ndims, dims[] = { (mwSize) numTimes, (mwSize) numPts, (mwSize) numPde };
      ndims = 3;
      typedef boost::const_multi_array_ref<double, 3> Matrix3;
      boost::array<Matrix3::index, 3> shape =
      { { numTimes, numPde, numPts } };
      Matrix3 sol3(pdeSol.getSolution().data(), shape, boost::fortran_storage_order());
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
      std::copy_n(pdeSol.getSolution().data(), numTimes*numPts*numPde, mxGetPr(sol));
    }

    int lhsIndex = 1;
    if (viewMesh > 1) {
      // return struct for results on a view mesh
      const char *fieldNames[] = { "x", "u", "uOde" };
      int numFields = 2;
      if (numOde)
        numFields++;
      mxArray *solStruc = mxCreateStructMatrix(1, 1, numFields, &fieldNames[0]);
      mxSetField(solStruc, 0, "x", MexInterface::toMxArray(pdeSol.getX().transpose()));
      mxSetField(solStruc, 0, "u", sol);
      if (numOde)
        mxSetField(solStruc, 0, "uOde", MexInterface::toMxArray(pdeSol.uOde));
      plhs[0] = solStruc;
    }
    else {
      plhs[0] = sol;
      // odes are included
      if (hasODE && nlhs > 1)
        plhs[lhsIndex++] = MexInterface::toMxArray(pdeSol.uOde);
    }

    // if we have events, more outputs are possible
    if (eventsFunc) {
      if (lhsIndex < nlhs) {
        // tsol
        plhs[lhsIndex++] = MexInterface::toMxArray(pdeSol.getOutputTimes());
      }
      if (lhsIndex < nlhs) {
        // sole
        plhs[lhsIndex++] = MexInterface::toMxArray(pdeSol.getEventsSolution());
      }
      if (lhsIndex < nlhs) {
        // te
        plhs[lhsIndex++] = MexInterface::toMxArray(pdeSol.getEventsTimes());
      }
      if (lhsIndex < nlhs) {
        // ie
        plhs[lhsIndex++] = MexInterface::toMxArray(pdeSol.getEventsIndex());
      }
    }

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

}