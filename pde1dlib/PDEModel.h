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
  typedef std::vector<int> DofList;
  void getDofIndicesForElem(int elemIndex,
    DofList &dofs) const;
  const PDEElement &element(int i) const { return elements[i]; }
  size_t numElements() const { return elements.size();  }
  size_t numEquations() const { return numEqns; }
  size_t numNodesFEEqns() const { return numNodesFEEqns_;  }
  template<class TG, class TM>
  void globalToMeshVec(const TG &gV, TM &mV);
private:
  const RealVector &origMesh;
  ShapeFunctionManager &sfm;
  std::vector<PDEElement> elements;
  std::vector<int> elementDofOffsets;
  size_t numNodesFEEqns_, numDofsPerNode;
  size_t numEqns; // total number of equations in the model
};

template<class TG, class TM>
void PDEModel::globalToMeshVec(const TG &gV, TM &mV) {
  size_t nnMesh = origMesh.rows();
  Eigen::Map<const RealMatrix> gM(gV.data(), numDofsPerNode, numNodesFEEqns_);
  Eigen::Map<RealMatrix> mM(mV.data(), numDofsPerNode, nnMesh);
  mM.col(0) = gM.col(0);
  int iOff = 0;
  for (int i = 0; i < numElements(); i++) {
    iOff += elements[i].numNodes() - 1;
    mM.col(i+1) = gM.col(iOff);
  }
}


