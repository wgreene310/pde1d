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



