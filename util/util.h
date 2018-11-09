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

#include <stdio.h>
#include <iostream>

#include <Eigen/Core>

#include <nvector/nvector_serial.h>

void print(N_Vector v, const char *title);

void pdePrintf(const char * format, ...);

Eigen::VectorXd linspace(double a, double b, int n);

template<class T>
void printMat(const T &a, const char *title, const char *format = "%16.9e,") {
  const size_t m = a.rows();
  const size_t n = a.cols();
  pdePrintf("%s(%zd,%zd)\n", title, m, n);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      pdePrintf(format, a(i, j));
    }
    pdePrintf("\n");
  }
  pdePrintf("\n");
}

template
<typename T, template<typename ELEM, typename ALLOC = std::allocator<ELEM>> class Container>
std::ostream& operator<< (std::ostream& out, const Container<T>& v)
{
  out << "{";
  //typename Container<T>::const_iterator beg = container.begin();
  const auto last = --v.end();
  for (auto it = v.begin(); it != v.end(); ++it) {
    out << *it;
    if (it != last)
      out << ", ";
  }
  out << "}";
  return out;
}