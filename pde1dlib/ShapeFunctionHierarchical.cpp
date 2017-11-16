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
 * ShapeFunctionHierarchical.cpp
 *
 *  Created on: Oct 30, 2016
 *      Author: bgreene
 */

#include <cmath>

#include "ShapeFunctionHierarchical.h"
#include "PDE1dException.h"

ShapeFunctionHierarchical::ShapeFunctionHierarchical(int degree) :
degree(degree)
{
  if (degree < 1 || degree > MAX_DEGREE)
    throw PDE1dException("pde1d:invalid_shape_degree",
  "ShapeFunctionHierarchical: invalid degree.");
  currentXi = -1;
}

void ShapeFunctionHierarchical::N(double xi, double *vals) const {
  vals[0] = (1 - xi) / 2.;
  vals[1] = (1 + xi) / 2.;
  bool xiOK = isXiOK(xi);
  if (degree >= 2) {
    if (!xiOK) {
      calcP(0, xi);
      calcP(1, xi);
    }
  }
  for(int i=2; i<=degree; i++) {
    if (!xiOK)
      calcP(i, xi);
    vals[i] = Ni(i+1);
  }
#if 0
  printf("P=");
  for (int i = 0; i <= degree ; i++)
    printf("%15.6e ", P[i]);
  printf("\n");
#endif
}

void ShapeFunctionHierarchical::dNdr(double xi, double *vals) const
{
  vals[0] = -.5;
  vals[1] = .5;
  bool xiOK = isXiOK(xi);
  if (degree >= 2) {
    if (!xiOK) {
      calcP(0, xi);
      calcP(1, xi);
    }
    dP[0] = 0;
    dP[1] = 1;
  }
  for (int i = 2; i <= degree; i++) {
    if (!xiOK)
      calcP(i, xi);
    vals[i] = dNi(i);
  }
#if 0
  printf("dP=");
  for (int i = 0; i <= degree; i++)
    printf("%15.6e ", dP[i]);
  printf("\n");
#endif
}

int ShapeFunctionHierarchical::numNodes() const {
  return degree + 1;
}



