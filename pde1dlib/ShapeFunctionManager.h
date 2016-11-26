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
    const Eigen::MatrixXd &N() const { return N_; }
    const Eigen::MatrixXd &dN() const { return dN_; }
    const Eigen::VectorXd &intRuleWts() const { return intRuleWts_;  }
    const ShapeFunction &getShapeFunction() const { return sf; }
  private:
    friend class ShapeFunctionManager;
    EvaluatedSF(int pOrder, int numIntPts);
    ShapeFunctionHierarchical sf;
    Eigen::MatrixXd N_, dN_;
    GausLegendreIntRule intRule;
    Eigen::VectorXd intRuleWts_;
  };
  const EvaluatedSF &getShapeFunction(int polyOrder);
private:
  typedef std::unordered_map<int, EvaluatedSF> SFMap;
  SFMap sfMap;
};

#endif /* PDE1DLIB_SHAPEFUNCTIONMANAGER_H_ */
