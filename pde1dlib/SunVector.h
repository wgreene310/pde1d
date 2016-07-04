// Copyright (C) 2016 William H. Greene
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

#ifndef SunVector_h
#define SunVector_h

#include <stddef.h>

#include <nvector/nvector_serial.h>

class SunVector {
public:
  SunVector(size_t n);
  ~SunVector();
  void setConstant(double c);
  double &operator[](int i) {
    return NV_DATA_S(nv)[i];
  }
  N_Vector operator()() {
    return nv;
  }
private:
  N_Vector nv;
};

#endif

