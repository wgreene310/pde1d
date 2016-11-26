#include "PDEElement.h"


PDEElement::PDEElement(int n1, int n2, 
  const ShapeFunctionManager::EvaluatedSF *esf) : esf(esf)
{
  nodeIndices[0] = n1;
  nodeIndices[1] = n2;
}


PDEElement::~PDEElement()
{
}
