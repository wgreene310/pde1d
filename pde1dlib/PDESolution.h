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
    return (int) time.size();
  }
  const RealVector &getX() const {
    return x;
  }
  const RealVector &getTime() const {
    return time;
  }
  void setSolutionVector(int timeStep, double time,
    const RealVector &u);
  const RealMatrix &getSolution() const {
    return u;
  }
  void print() const;
//private:
  RealMatrix uOde;
private:
  const RealVector &initialX;
  const RealVector &time;
  RealVector x;
  int numViewElemsPerElem;
  size_t numNodesViewMesh;
  const PDE1dDefn &pde;
  const PDEModel &model;
  RealMatrix u;
  // temporary work storage;
  RealMatrix u2Tmp;
};

#endif /* PDE1DLIB_PDESOLUTION_H_ */
