#pragma once

#include "ShapeFunctionManager.h"

class PDEElement
{
public:
  PDEElement(int n1 = 0, int n2 = 0, const ShapeFunctionManager::EvaluatedSF *esf= 0);
  ~PDEElement();
  const ShapeFunctionManager::EvaluatedSF &getSF() const { return *esf; }
private:
  const ShapeFunctionManager::EvaluatedSF *esf;
  int nodeIndices[2];
};

