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

#include "PDE1dDefn.h"
#include "PDE1dImpl.h"
#include "PDE1dOptions.h"


PDE1dDefn::PDE1dDefn()
{
}

PDESolution pde1d(PDE1dDefn &pde)
{
  PDE1dOptions options;
  PDE1dImpl pdeImpl(pde, options);
  PDESolution sol;
  int err = pdeImpl.solveTransient(sol);
  return sol;
}