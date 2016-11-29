#pragma once

#include "ShapeFunctionManager.h"

class PDEElement
{
public:
  PDEElement(int n1, int n2, const ShapeFunctionManager::EvaluatedSF *esf);
  PDEElement(); // required for std::vector
  ~PDEElement();
  const ShapeFunctionManager::EvaluatedSF &getSF() const { return *esf; }
  int numNodes() const { return numNodes_;  }
private:
  const ShapeFunctionManager::EvaluatedSF *esf;
  int nodeIndices[2];
  int numNodes_;
};

