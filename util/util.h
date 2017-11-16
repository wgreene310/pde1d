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

#include <Eigen/Core>

#include <nvector/nvector_serial.h>

void print(N_Vector v, const char *title);

Eigen::VectorXd linspace(double a, double b, int n);

template<class T>
void printMat(const T &a, const char *title, const char *format = "%16.9e,") {
  const size_t m = a.rows();
  const size_t n = a.cols();
  printf("%s(%zd,%zd)\n", title, m, n);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      printf(format, a(i, j));
    }
    printf("\n");
  }
  printf("\n");
}