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
 * PDESolution.h
 *
 *  Created on: Apr 30, 2017
 *      Author: bgreene
 */

#ifndef PDE1DLIB_PDESOLUTION_H_
#define PDE1DLIB_PDESOLUTION_H_

#include <MatrixTypes.h>

class PDE1dDefn;
class PDEModel;

class PDESolution {
public:
  PDESolution(const PDE1dDefn &pde, const PDEModel &model,
    int numViewElemsPerElem);
  int numSpatialPoints() const {
    return (int) x.size();
  }
  int numTimePoints() const {
    return (int) outTimes.size();
  }
  const RealVector &getX() const {
    return x;
  }
  const RealVector &getOutputTimes() const {
    return outTimes;
  }
  void setSolutionVector(int timeStep, double time,
    const RealVector &u);
  const RealMatrix &getSolution() const {
    return u;
  }
  void print() const;
  void close();
  void setEventsSolution(int timeStep, double time,
    const RealVector &u, const IntVector &eventsFound);
  const RealMatrix getEventsSolution() const {
    return eventsSolution;
  }
  const RealVector getEventsTimes() const {
    return eventsTimes;
  }
  const IntVector getEventsIndex() const {
    return eventsIndex;
  }
//private:
  RealMatrix uOde;
private:
  int numEvents, numEventsSet;
  int numTimesSet;
  const RealVector &initialX;
  const RealVector &time;
  RealVector x;
  int numViewElemsPerElem;
  size_t numNodesViewMesh;
  const PDE1dDefn &pde;
  const PDEModel &model;
  RealMatrix u;
  RealVector outTimes;
  RealMatrix u2Tmp; // temporary work storage;
  RealMatrix eventsSolution;
  RealVector eventsTimes;
  IntVector eventsIndex;
};

#endif /* PDE1DLIB_PDESOLUTION_H_ */
