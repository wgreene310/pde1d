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

#include <iostream>
#include <stdio.h>
#include <cstdarg>

#include "util.h"

void print(N_Vector v, const char *title) {
  double *d = NV_DATA_S(v);
  int len = NV_LENGTH_S(v);
  int count = 0;
  pdePrintf("Vector: %s(%d)\n", title, len);
  for (int i = 0; i < len; i++) {
    pdePrintf(" %14.9e", d[i]);
    if (++count == 6) {
      pdePrintf("\n");
      count = 0;
    }
  }
  if (count)
    pdePrintf("\n");
}

#if 0
void pdePrintf(const char * format, ...)
{
  va_list ap;
  va_start(ap, format);
  char msg[4096];
  vsprintf(msg, format, ap);
  va_end(ap);
  // matlab mex connects cout to console but not printf
  //std::cout << msg;
  printf(msg);
}
#endif

Eigen::VectorXd linspace(double start, double end, int n)
{
  Eigen::VectorXd v(n);
  int nm1 = n - 1;
  double dx = (end - start) / nm1;
  double x = start;
  for (int i = 0; i < n; i++) {
    v[i] = x;
    x += dx;
  }
  return v;
}