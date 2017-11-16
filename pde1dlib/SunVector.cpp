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

#include "SunVector.h"
#include "PDE1dException.h"


SunVector::SunVector(size_t n) : 
_SundialsVector_(N_VNew_Serial((int) n)),
Eigen::Map<Eigen::VectorXd>(NV_DATA_S(nv), n)
{
  if (!nv) {
    char msg[256];
    sprintf(msg, "Sundials error: unable to allocate serial vector of "
      "length %zu", n);
    throw PDE1dException("pde1d:sundials_mem_alloc", msg);
  }
  isExternal = false;
}

SunVector::SunVector(const SunVector &rhs) : 
_SundialsVector_(N_VNew_Serial((int) rhs.size())), 
Eigen::Map<Eigen::VectorXd>(NV_DATA_S(nv), rhs.size())
{
  if (!nv) {
    char msg[256];
    sprintf(msg, "Sundials error: unable to allocate serial vector of "
      "length %zu", size());
    throw PDE1dException("pde1d:sundials_mem_alloc", msg);
  }
  Eigen::Map<Eigen::VectorXd>::operator=(rhs);
  isExternal = false;
}

SunVector::SunVector(N_Vector nv) : _SundialsVector_(nv),
Eigen::Map<Eigen::VectorXd>(NV_DATA_S(nv), NV_LENGTH_S(nv))
{
  isExternal = true;
}

SunVector::~SunVector()
{
  if (nv && !isExternal)
    N_VDestroy_Serial(nv);
}

