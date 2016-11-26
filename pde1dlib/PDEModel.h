#pragma once
#include <vector>

#include <MatrixTypes.h>
#include "PDEElement.h"

class ShapeFunctionManager;

class PDEModel
{
public:
  PDEModel(const RealVector &mesh, int polyOrder, size_t numDofsPerNode,
    ShapeFunctionManager &sfm);
  ~PDEModel();
  void getDofIndicesForElem(int elemIndex,
    std::vector<int> &dofs) const;
  const PDEElement &element(int i) const { return elements[i]; }
  const size_t numElements() const { return elements.size();  }
private:
  const RealVector &origMesh;
  ShapeFunctionManager &sfm;
  std::vector<PDEElement> elements;
  std::vector<int> elementDofOffsets;
  size_t numDofsPerNode;
};

