#include "PDEModel.h"
#include "ShapeFunctionManager.h"


PDEModel::PDEModel(const RealVector &mesh, int polyOrder, size_t numDofsPerNode,
  ShapeFunctionManager &sfm) :
origMesh(mesh), numDofsPerNode(numDofsPerNode),
sfm(sfm)
{
  size_t numElems = mesh.size()-1;
  elements.resize(numElems);
  elementDofOffsets.resize(numElems);
  const ShapeFunctionManager::EvaluatedSF &esf =
    sfm.getShapeFunction(polyOrder);
  int dof = 0;
  for (int i = 0; i < numElems; i++) {
    elements[i] = PDEElement(i, i + 1, &esf);
    elementDofOffsets[i] = dof;
    dof += polyOrder;
  }
}


PDEModel::~PDEModel()
{
}

void PDEModel::getDofIndicesForElem(int elemIndex,
  std::vector<int> &dofs) const
{
  const PDEElement &elem = elements[elemIndex];
  const ShapeFunctionManager::EvaluatedSF &esf = elem.getSF();
  int nn = esf.getShapeFunction().numNodes();
  dofs.resize(nn);
  int d = elementDofOffsets[elemIndex];
  dofs[0] = d;
  dofs[1] = d + nn - 1;
  for (int i = 2; i < nn; i++)
    dofs[i] = ++d;
}