
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
