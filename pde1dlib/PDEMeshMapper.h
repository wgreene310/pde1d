// Copyright (C) 2016-2017 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

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

