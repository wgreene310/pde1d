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

#ifndef PDE1dMexInt_h
#define PDE1dMexInt_h

#include <algorithm>

#include <mex.h>

#include <PDE1dDefn.h>

#ifdef _MSC_VER
// octave incorrectly defines mwSize
#pragma warning( push )
#pragma warning( disable : 4267)
#endif

class PDE1dMexInt : public PDE1dDefn
{
public:
  PDE1dMexInt(int m, const mxArray *pdefun, const mxArray *icfun,
    const mxArray *bcfun,
    const mxArray *xmesh, const mxArray *tspan);
  ~PDE1dMexInt();
  //virtual int getNumEquations();
  virtual int getCoordSystem() const {
    return mCoord;
  }
  size_t numNodes() { return mesh.size(); }
  virtual void evalIC(double x, RealVector &ic);
  virtual void evalODEIC(RealVector &ic);
  virtual void evalBC(double xl, const RealVector &ul,
    double xr, const RealVector &ur, double t, 
    const RealVector &v, const RealVector &vDot, BC &bc);
  virtual void evalPDE(double x, double t,
    const RealVector &u, const RealVector &DuDx, 
    const RealVector &v, const RealVector &vDot, PDE &pde);
  virtual bool hasVectorPDEEval() const { return true; }
  virtual void evalPDE(const RealVector &x, double t,
    const RealMatrix &u, const RealMatrix &DuDx, 
    const RealVector &v, const RealVector &vDot, PDECoeff &pde);
  virtual const RealVector &getMesh() const;
  virtual const RealVector &getODEMesh();
  virtual const RealVector &getTimeSpan() const;
  void setODEDefn(const mxArray *odeFun, const mxArray *icFun,
    const mxArray *odeMesh);
  virtual int getNumODE() const { return numODE; }
  virtual int getNumPDE() const { return numPDE; }
  virtual void evalODE(double t, const RealVector &v,
    const RealVector &vdot,
    const RealMatrix &u, const RealMatrix &DuDx, const RealMatrix &R,
    const RealMatrix &odeDuDt, const RealMatrix &odeDuDxDt,
    RealVector &f);
  void setEventsFunction(const mxArray *eventsFun);
  virtual int getNumEvents() const { return numEvents; }
  virtual void evalEvents(double t, const RealMatrix &u,
    RealVector &eventsVal, RealVector &eventsIsTerminal,
    RealVector &eventsDirection);
private:
  static void setScalar(double x, mxArray *a) {
    const int vr = 1, vc = 1;
    if (vr != mxGetM(a) || vc != mxGetN(a)) {
      double *pr = (double*)mxRealloc(mxGetPr(a), sizeof(double));
      mxSetM(a, vr);
      mxSetN(a, vc);
      mxSetPr(a, pr);
    }
    double *p = mxGetPr(a); p[0] = x;
  }
  template<class T>
  static void setMxImpl(const T &ea, mxArray *a) {
    size_t vr = ea.rows(), vc = ea.cols(), s = ea.size();
    if (vr != mxGetM(a) || vc != mxGetN(a)) {
      double *pr = (double*)mxRealloc(mxGetPr(a), s*sizeof(double));
      mxSetM(a, vr);
      mxSetN(a, vc);
      mxSetPr(a, pr);
    }
    std::copy_n(ea.data(), s, mxGetPr(a));
  }
  template<class T>
  static void setVector(const T &v, mxArray *a) {
    size_t s = v.size();
    if (s != mxGetM(a) || mxGetN(a) != 1) {
      double *pr = (double*)mxRealloc(mxGetPr(a), s * sizeof(double));
      mxSetM(a, s);
      mxSetN(a, 1);
      mxSetPr(a, pr);
    }
    std::copy_n(v.data(), s, mxGetPr(a));
  }
  static void setMatrix(const RealMatrix &v, mxArray *a) {
    setMxImpl(v, a);
  }
  void setNumPde();
  void setNumOde();
  void setNumEvents();
  void callMatlab(const mxArray *inArgs[], int nargin);
  void callMatlab(const mxArray *inArgs[], int nargin,
    RealVector *outArgs[], int nargout);
  void callMatlab(const mxArray *inArgs[], int nargin,
    RealMatrix *outArgs[], int nargout);
  void callMatlab(const mxArray *inArgs[], int nargin, int nargout);
  static std::string getFuncNameFromHandle(const mxArray *funcHandle);
  void checkMxType(const mxArray *a, int argIndex, const mxArray *funcHandle);
  static void destroy(mxArray *a) {
    // not sure if mex does this test?
    if (a)
      mxDestroyArray(a);
  }
  int mCoord, numPDE, numODE, numEvents;
  RealVector mesh, tSpan, odeMeshVec;
  const mxArray *pdefun, *icfun, *bcfun, *xmesh, *tspan;
  const mxArray *odefun, *odeIcFun, *odeMesh;

  mxArray *mxX1, *mxX2, *mxT; // scalar x and t
  mxArray *mxVec1, *mxVec2, *mxMat1, *mxMat2; // input to pde function
  mxArray *mxV, *mxVDot, *mxOdeU, *mxOdeDuDx, *mxOdeR, 
    *mxOdeDuDt, *mxOdeDuDxDt;
  static const int maxMatlabRetArgs = 4;
  mxArray *matOutArgs[maxMatlabRetArgs];

  const mxArray *eventsFun;
  mxArray *mxM;
  mxArray *mxEventsU; // solution at mesh points
};

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif

