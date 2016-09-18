#pragma once

#include <MatrixTypes.h>

class ShapeFunction;

class PDEMeshMapper {
public:
  PDEMeshMapper(const RealVector &srcMesh, const ShapeFunction &sf,
    const RealVector &destMesh);
  ~PDEMeshMapper();
  void mapFunction(const RealMatrix &srcU, RealMatrix &destU);
private:
  const RealVector &srcMesh, &destMesh;
  const ShapeFunction &sf;
  RealVector destMeshParamVals;
  Eigen::VectorXi destMeshElemIndex;
};

