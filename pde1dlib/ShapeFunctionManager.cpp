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

/*
 * ShapeFunctionManager.cpp
 *
 *  Created on: Nov 12, 2016
 *      Author: bgreene
 */

#include "ShapeFunctionManager.h"

ShapeFunctionManager::EvaluatedSF::EvaluatedSF(int pOrder, int numIntPts) :
intRule(GausLegendreIntRule(numIntPts)), 
sf(ShapeFunctionHierarchical(pOrder))
{
}


ShapeFunctionManager::ShapeFunctionManager()
{

}

const ShapeFunctionManager::EvaluatedSF &ShapeFunctionManager::getShapeFunction(int polyOrder)
{
  int numIntPts = GausLegendreIntRule::getNumPtsForPolyOrder(2 * polyOrder);
  SFMap::iterator it = sfMap.find(polyOrder);
  if (it == sfMap.end()) {
    EvaluatedSF esf(polyOrder, numIntPts);
    int ndof = polyOrder + 1;
    esf.N_.resize(ndof, numIntPts);
    esf.dN_.resize(ndof, numIntPts);
    esf.intRuleWts_.resize(numIntPts);
    for (int i = 0; i < numIntPts; i++) {
      double pt;
      esf.intRule.getPoint(i, pt, esf.intRuleWts_[i]);
      esf.sf.N(pt, esf.N_.col(i).data());
      esf.sf.dNdr(pt, esf.dN_.col(i).data());
    }
    auto iit = sfMap.insert(std::make_pair(polyOrder, esf));
    it = iit.first;
  }
  return it->second;
}