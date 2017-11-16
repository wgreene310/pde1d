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
using std::cout;
using std::endl;

#include "PDEElement.h"

PDEElement::PDEElement(int n1, int n2, 
  const ShapeFunctionManager::EvaluatedSF *esf) : esf(esf)
{
  nodeIndices[0] = n1;
  nodeIndices[1] = n2;
  numNodes_ = esf->getShapeFunction().numNodes();
}

PDEElement::PDEElement()
{
  esf = 0;
  numNodes_ = 0;
  nodeIndices[0] = 0;
  nodeIndices[1] = 0;
}


PDEElement::~PDEElement()
{
}
