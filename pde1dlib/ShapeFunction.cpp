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

#include "ShapeFunction.h"

void ShapeFunction2::N(double r, double *func) const {
  func[0] = .5*(1 - r);
  func[1] = .5*(1 + r);
}

void ShapeFunction2::dNdr(double r, double *df) const {
  df[0] = -.5;
  df[1] = .5;
}

void ShapeFunction3::N(double r, double *func) const {
  func[0] = -r*(1 - r) / 2.;
  func[1] = r*(1 + r) / 2.;
  func[2] = 1 - r*r;
}

void ShapeFunction3::dNdr(double r, double *df) const {
  df[0] = (-1 + 2 * r) / 2.;
  df[1] = (1 + 2 * r) / 2.;
  df[2] = -2 * r;
}