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

/*
 * PDESolution.cpp
 *
 *  Created on: Apr 30, 2017
 *      Author: bgreene
 */
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

#include "PDESolution.h"
#include "PDE1dDefn.h"
#include "PDEModel.h"

PDESolution::PDESolution(const PDE1dDefn &pde, const PDEModel &model,
  int numViewElemsPerElem) :
  initialX(pde.getMesh()),
  time(pde.getTimeSpan()), numViewElemsPerElem(numViewElemsPerElem),
  pde(pde), model(model)
{
  numTimesSet = 0;
  numEventsSet = 0;
  size_t maxNumTimes = time.size();
  outTimes.resize(maxNumTimes);
  int numPDE = pde.getNumPDE();
  if (numViewElemsPerElem == 1) {
    x = initialX;
    size_t nnOrig = x.size();
    u.resize(maxNumTimes, nnOrig*numPDE);
    const size_t nnfee = model.numNodesFEEqns();
    u2Tmp.resize(1, nnOrig*numPDE);
  }
  else {
    size_t ne = initialX.size() - 1;
    numNodesViewMesh = numViewElemsPerElem*ne + 1;
    x.resize(numNodesViewMesh);
    u.resize(maxNumTimes, numNodesViewMesh*numPDE);
    u2Tmp.resize(numPDE, numNodesViewMesh);
    x(0) = initialX(0);
    double xij = x(0);
    size_t ij = 1;
    for (int i = 0; i < ne; i++) {
      double L = initialX(i + 1) - initialX(i);
      double dx = L / (double)numViewElemsPerElem;
      for (int j = 0; j < numViewElemsPerElem; j++) {
        xij += dx;
        x(ij++) = xij;
      }
    }
  }
  numEvents = pde.getNumEvents();
  if (numEvents) {
    eventsSolution.resize(numEvents, u.cols());
    eventsTimes.resize(numEvents);
    eventsIndex.resize(numEvents);
  }
}

void PDESolution::setSolutionVector(int timeStep, double time,
  const RealVector &uSol)
{
  outTimes(numTimesSet++) = time;
  int numDepVars = pde.getNumPDE();
  const size_t nnfee = model.numNodesFEEqns();
  if (numViewElemsPerElem == 1) {
    model.globalToMeshVec(uSol.topRows(nnfee), u2Tmp);
    u.row(timeStep) = u2Tmp;
#if 0
    int numDepVars = pde.getNumPDE();
    const size_t nnfee = model.numNodesFEEqns();
    typedef Eigen::Map<const Eigen::MatrixXd> ConstMapMat;
    ConstMapMat u2Sol(uSol.data(), numDepVars, nnfee);
    cout << "u2Sol=" << u2Sol << endl;
#endif
  }
  else {
    typedef Eigen::Map<const Eigen::MatrixXd> ConstMapMat;
    ConstMapMat u2Sol(uSol.data(), numDepVars, nnfee);
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
    //cout << "rvals=" << rVals << endl;
    RealMatrix N, u2eView;
    size_t ne = model.numElements();
    int lastNumElemNodes = 0;
    size_t iViewNode = 0;
    for (int e = 0; e < ne; e++) {
      const PDEElement &elem = model.element(e);
      model.getDofIndicesForElem(e, eDofs);
      PDEModel::globalToElemVec(eDofs, u2Sol, u2e);
      const ShapeFunction &sf = elem.getSF().getShapeFunction();
      if (sf.numNodes() != lastNumElemNodes) {
        lastNumElemNodes = sf.numNodes();
        N.resize(lastNumElemNodes, numR);
        for (int i = 0; i < numR; i++)
          sf.N(rVals(i), N.col(i).data());
        //cout << "N\n" << N << endl;
      }
      u2eView = u2e*N;
      u2Tmp.col(iViewNode) = u2e.col(0);
      u2Tmp.block(0, iViewNode + 1, numDepVars, numR) = u2eView;
      iViewNode += numR + 1;
    }
    // finish the last view node
    u2Tmp.col(iViewNode) = u2e.col(1);
    MapVec u2ViewVector(u2Tmp.data(), u2Tmp.size());
    u.row(timeStep) = u2ViewVector;
  }
}

void PDESolution::print() const
{
  cout << "Solution, u\n" << u << endl;
}

void PDESolution::close()
{
  outTimes.conservativeResize(numTimesSet);
  u.conservativeResize(numTimesSet, u.cols());
  if (numEventsSet < numEvents) {
    eventsSolution.conservativeResize(numEventsSet, u.cols());
    eventsTimes.conservativeResize(numEventsSet);
    eventsIndex.conservativeResize(numEventsSet);
  }
}

void PDESolution::setEventsSolution(int timeStep, double time,
  const RealVector &u, const IntVector &eventsFound)
{
  eventsSolution.row(numEventsSet) = u;
  eventsTimes(numEventsSet) = time;
  eventsIndex(numEventsSet) = numEventsSet + 1;
  ++numEventsSet;
}