/*
 * ShapeFunctionManager.cpp
 *
 *  Created on: Nov 12, 2016
 *      Author: bgreene
 */

#include "ShapeFunctionManager.h"

ShapeFunctionManager::ShapeFunctionManager()
{

}

const ShapeFunctionManager::EvaluatedSF &ShapeFunctionManager::getShapeFunction(int polyOrder, 
  int numIntPts)
{
  IIP p(polyOrder, numIntPts);
  SFMap::iterator it = sfMap.find(p);
  if (it == sfMap.end()) {
    EvaluatedSF esf(numIntPts);
    int ndof = polyOrder + 1;
    esf.N_.resize(ndof, numIntPts);
    esf.dN_.resize(ndof, numIntPts);
    esf.intRuleWts_.resize(numIntPts);
    ShapeFunctionHierarchical sf(polyOrder);
    for (int i = 0; i < numIntPts; i++) {
      double pt;
      esf.intRule.getPoint(i, pt, esf.intRuleWts_[i]);
      sf.N(pt, esf.N_.col(i).data());
      sf.dNdr(pt, esf.dN_.col(i).data());
    }
    auto iit = sfMap.insert(std::make_pair(p, esf));
    it = iit.first;
  }
  return it->second;
}