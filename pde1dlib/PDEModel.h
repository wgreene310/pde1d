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
  void globalToMeshVec(const TG &gV, TM &mV) const;
  static void globalToElemVec(const DofList &eDofs, const RealMatrix &ug,
    RealMatrix &ue);
  template<class TE, class TG>
  static void elemToGlobalVec(const PDEModel::DofList &eDofs, const TE &ue, TG &ug);
  template<class TE, class TG>
  static void assembleElemVec(const PDEModel::DofList &eDofs, size_t numNds, size_t numDofPerNode,
    const TE &eVec, TG &gVec);
private:
  const RealVector &origMesh;
  ShapeFunctionManager &sfm;
  std::vector<PDEElement> elements;
  std::vector<int> elementDofOffsets;
  size_t numNodesFEEqns_, numDofsPerNode;
  size_t numEqns; // total number of equations in the model
};

template<class TG, class TM>
void PDEModel::globalToMeshVec(const TG &gV, TM &mV) const {
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

template<class TE, class TG>
void PDEModel::elemToGlobalVec(const DofList &eDofs, const TE &ue, TG &ug) {
  size_t n = eDofs.size();
  for (size_t i = 0; i < n; i++)
    ug.col(eDofs[i]) = ue.col(i);
}

template<class TE, class TG>
void PDEModel::assembleElemVec(const DofList &eDofs, size_t numNds, 
  size_t numDofPerNode, const TE &eVec, TG &gVec) {
  size_t n = eDofs.size();
  Eigen::Map<const RealMatrix> eMat(eVec.data(), numDofPerNode, n);
  Eigen::Map<RealMatrix> gMat(&gVec[0], numDofPerNode, numNds);
  for (size_t i = 0; i < n; i++)
    gMat.col(eDofs[i]) += eMat.col(i);
}


