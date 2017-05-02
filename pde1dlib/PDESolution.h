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
