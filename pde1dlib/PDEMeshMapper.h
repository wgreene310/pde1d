#pragma once

#include <MatrixTypes.h>

class PDEModel;

class PDEMeshMapper {
public:
  PDEMeshMapper(const RealVector &srcMesh, const PDEModel &model,
    const RealVector &destMesh);
  void mapFunction(const RealMatrix &srcU, RealMatrix &destU);
  void mapFunctionDer(const RealMatrix &srcU, RealMatrix &destU);
private:
  void mapFunctionImpl(const RealMatrix &srcU, RealMatrix &destU, bool calcDeriv);
  const RealVector &srcMesh, &destMesh;
  const PDEModel &model;
  RealVector destMeshParamVals;
  Eigen::VectorXi destMeshElemIndex;
};

