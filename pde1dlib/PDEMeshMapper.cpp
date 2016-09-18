#include "PDEMeshMapper.h"


PDEMeshMapper::PDEMeshMapper(const RealVector &srcMesh, const ShapeFunction &sf,
  const RealVector &destMesh) :
  srcMesh(srcMesh), destMesh(destMesh), sf(sf)
{
  const int numDest = destMesh.size();
  destMeshParamVals.resize(numDest);
  destMeshElemIndex.resize(numDest);
}


PDEMeshMapper::~PDEMeshMapper()
{
}

void PDEMeshMapper::mapFunction(const RealMatrix &srcU, RealMatrix &destU)
{

}
