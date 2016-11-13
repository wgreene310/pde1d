/*
 * ShapeFunctionManager.h
 *
 *  Created on: Nov 12, 2016
 *      Author: bgreene
 */

#ifndef PDE1DLIB_SHAPEFUNCTIONMANAGER_H_
#define PDE1DLIB_SHAPEFUNCTIONMANAGER_H_

#include <unordered_map>

#include <Eigen/Core>

#include "ShapeFunctionHierarchical.h"
#include "GausLegendreIntRule.h"

class ShapeFunctionManager {
public:
  ShapeFunctionManager();
  class EvaluatedSF {
  public:
    EvaluatedSF(int numIntPts) : 
      intRule(GausLegendreIntRule(numIntPts)) {}
    const Eigen::MatrixXd &N() const { return N_; }
    const Eigen::MatrixXd &dN() const { return dN_; }
    const Eigen::VectorXd &intRuleWts() const { return intRuleWts_;  }
  private:
    friend class ShapeFunctionManager;
    Eigen::MatrixXd N_, dN_;
    GausLegendreIntRule intRule;
    Eigen::VectorXd intRuleWts_;
  };
  const EvaluatedSF &getShapeFunction(int polyOrder, int numIntpts);
private:
  typedef std::pair<int, int> IIP;
  struct HashIntPair {
    size_t operator()(const IIP &p) const {
      return std::hash < int > {}(10 * p.first + p.second);
    }
  };
  typedef std::unordered_map<IIP, EvaluatedSF, HashIntPair> SFMap;
  SFMap sfMap;
};

#endif /* PDE1DLIB_SHAPEFUNCTIONMANAGER_H_ */
