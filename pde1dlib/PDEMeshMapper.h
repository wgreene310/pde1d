#pragma once

#include <MatrixTypes.h>

class ShapeFunction;

class PDEMeshMapper {
public:
  PDEMeshMapper(const RealVector &srcMesh, const ShapeFunction &sf,
    const RealVector &destMesh);
  void mapFunction(const RealMatrix &srcU, RealMatrix &destU);
  void mapFunctionDer(const RealMatrix &srcU, RealMatrix &destU);
private:
  void mapFunctionImpl(const RealMatrix &srcU, RealMatrix &destU, bool calcDeriv);
  const RealVector &srcMesh, &destMesh;
  const ShapeFunction &sf;
  RealVector destMeshParamVals;
  Eigen::VectorXi destMeshElemIndex;
};

