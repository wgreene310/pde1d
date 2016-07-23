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

#include <stdio.h>

#include "SunVector.h"
#include "PDE1dException.h"


SunVector::SunVector(size_t n)
{
  nv = N_VNew_Serial(n);
  if (!nv) {
    char msg[256];
    sprintf(msg, "Sundials error: unable to allocate serial vector of "
      "length %d", n);
    throw PDE1dException("pde1d:sundials_mem_alloc", msg);
  }
}


SunVector::~SunVector()
{
  if (nv) 
    N_VDestroy_Serial(nv);
}

void SunVector::setConstant(double c)
{
  N_VConst(c, nv);
}
