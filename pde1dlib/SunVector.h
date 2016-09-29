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
#include <assert.h>

#include <Eigen/Core>

#include <nvector/nvector_serial.h>

struct _SundialsVector_ {
  _SundialsVector_(N_Vector nv) : nv(nv) {}
  N_Vector nv;
};


class SunVector : private _SundialsVector_ ,
  public Eigen::Map<Eigen::VectorXd>
{
public:
  explicit SunVector(size_t n);
  explicit SunVector(N_Vector nv);
  SunVector(const SunVector &rhs);
  SunVector &operator=(const SunVector &rhs) {
    if (this == &rhs) return *this;
    assert(rows() == rhs.rows());
    Eigen::Map<Eigen::VectorXd>::operator=(rhs);
    return *this;
  }
  SunVector &operator=(const Eigen::Ref<const Eigen::VectorXd> &rhs)
  {
    assert(rows() == rhs.rows());
    Eigen::Map<Eigen::VectorXd>::operator=(rhs);
    return *this;
  }
  ~SunVector();
  N_Vector getNV() {
    return nv;
  }
private:
  bool isExternal;
};

#endif

