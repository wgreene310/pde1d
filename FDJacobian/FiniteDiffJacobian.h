// Copyright (C) 2016 William H. Greene
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

#include <nvector/nvector_serial.h>
#include <ida/ida.h>
#include <ida/ida_klu.h>

#include <Eigen/SparseCore>

struct SunSparseMat : public _SlsMat {
  SunSparseMat(Eigen::SparseMatrix<double> &eMat) {
    M = static_cast<int>(eMat.rows());
    N = static_cast<int>(eMat.cols());
    NNZ = static_cast<int>(eMat.nonZeros());
    data = eMat.valuePtr();
    rowvals = eMat.innerIndexPtr();
    colptrs = eMat.outerIndexPtr();
  }
};

class FiniteDiffJacobian {
public:
  typedef Eigen::SparseMatrix<double> SparseMat;
  FiniteDiffJacobian(SparseMat &jacPattern);
  ~FiniteDiffJacobian();
  void calcJacobian(double tres, double alpha, double beta,
    N_Vector uu, N_Vector up, N_Vector r,
    IDAResFn rf, void *userData, SparseMat &Jac, bool useCD=false);
  void calcJacobian(double tres, double alpha, double beta,
    N_Vector uu, N_Vector up, N_Vector r,
    IDAResFn rf, void *userData, SlsMat Jac);
  void calcJacobianCD(double tres, double alpha, double beta,
    N_Vector uu, N_Vector up, N_Vector r,
    IDAResFn rf, void *userData, SlsMat Jac);
private:
  int neq, nnz;
  Eigen::VectorXi indrow, jpntr, ngrp;
  int maxgrp, mingrp;
};

