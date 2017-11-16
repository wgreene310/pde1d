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
 * ShapeFunctionHierarchical.h
 *
 *  Created on: Oct 30, 2016
 *      Author: bgreene
 */

#ifndef PDE1DLIB_SHAPEFUNCTIONHIERARCHICAL_H_
#define PDE1DLIB_SHAPEFUNCTIONHIERARCHICAL_H_

#include <boost/math/special_functions/legendre.hpp>

#include "ShapeFunction.h"

class ShapeFunctionHierarchical: public ShapeFunction {
public:
  ShapeFunctionHierarchical(int degree);
  virtual void N(double r, double *vals) const;
  virtual void dNdr(double r, double *vals) const;
  virtual int numNodes() const;
  static const int MAX_DEGREE = 9;
private:
  int degree;
  mutable double P[MAX_DEGREE+1], dP[MAX_DEGREE+1];
  mutable double currentXi;
  inline double Ni(int i) const {
    return (P[i-1] - P[i - 3]) / sqrt(2. * (2. * (i - 1) - 1.));
  }
  inline double dPi(int i) const {
    return (2. * i + 1)*P[i] + dP[i - 1];
  }
  inline double dNi(int i) const {
    dP[i] = dPi(i-1);
    return (dP[i] - dP[i - 2]) / sqrt(2. * (2. * i - 1.));
  }
  bool isXiOK(double xi) const {
    if (currentXi == xi)
      return true;
    currentXi = xi;
    return false;
  }
  void calcP(int i, double xi) const {
    P[i] = boost::math::legendre_p(i, xi);
  }
};

#endif /* PDE1DLIB_SHAPEFUNCTIONHIERARCHICAL_H_ */
