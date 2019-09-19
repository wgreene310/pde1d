// Copyright (C) 2016-2019 William H. Greene
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
#include <stdio.h>
#include <algorithm>

#include "FiniteDiffJacobian.h"
#include "FDJacobian.h"
#include "SunVector.h"

namespace {
  const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
  inline double stepLen(double vi) {
    return sqrtEps*std::max(std::abs(vi), 1.0);
  }
}

FiniteDiffJacobian::FiniteDiffJacobian(SparseMat &jacPattern)
{
  neq = static_cast<int>(jacPattern.rows());
  nnz = static_cast<int>(jacPattern.nonZeros());
  indrow.resize(nnz);
  jpntr.resize(nnz);
  ngrp.resize(neq);
  Eigen::VectorXi indcol(nnz);
  int ii = 0;
  for (int k = 0; k < jacPattern.outerSize(); ++k)
    for (SparseMat::InnerIterator it(jacPattern, k); it; ++it) {
    indrow[ii] = static_cast<int>(it.row() + 1);   // row index
    indcol[ii] = static_cast<int>(it.col() + 1);
    ii++;
    }
  const int liwa = 6 * neq;
  Eigen::VectorXi ipntr(neq + 1), iwa(liwa);
  int info;
  dsm_(&neq, &neq, &nnz, indrow.data(), indcol.data(), ngrp.data(),
    &maxgrp, &mingrp, &info, ipntr.data(), jpntr.data(), iwa.data(), &liwa);
  //printf("info=%d, maxgrp=%d, mingrp=%d\n", info, maxgrp, mingrp);
}


FiniteDiffJacobian::~FiniteDiffJacobian()
{
}

void FiniteDiffJacobian::calcJacobian(double tres, double alpha,
  double beta, N_Vector uu, N_Vector up, N_Vector r,
  IDAResFn rf, void *userData, SparseMat &jac, bool useCD)
{
  jac.resize(neq, neq);
  jac.reserve(nnz);
  if (useCD)
    calcJacobianCD(tres, alpha, beta, uu, up, r, rf, userData,
      jac.valuePtr(), jac.outerIndexPtr(), jac.innerIndexPtr());
  else
    calcJacobian(tres, alpha, beta, uu, up, r, rf, userData,
      jac.valuePtr(), jac.outerIndexPtr(), jac.innerIndexPtr());
  jac.resizeNonZeros(nnz);
}

void FiniteDiffJacobian::calcJacobian(double tres, double alpha, double beta,
  N_Vector uu, N_Vector up, N_Vector r,
  IDAResFn rf, void *userData, SparseMap &jac, bool useCD) {
  if (useCD)
    calcJacobianCD(tres, alpha, beta, uu, up, r, rf, userData,
      jac.valuePtr(), jac.outerIndexPtr(), jac.innerIndexPtr());
  else
    calcJacobian(tres, alpha, beta, uu, up, r, rf, userData,
      jac.valuePtr(), jac.outerIndexPtr(), jac.innerIndexPtr());
}

void FiniteDiffJacobian::calcJacobian(double tres, double alpha,
  double beta, N_Vector uu, N_Vector up, N_Vector r,
  IDAResFn rFunc, void *userData, double *jacVals,
  int *jacColPtrs, int *jacRowIndices)
{
  typedef Eigen::Map<Eigen::VectorXd> MapVec;
  MapVec u(NV_DATA_S(uu), neq);
  MapVec uDot(NV_DATA_S(up), neq);
  MapVec resvec(NV_DATA_S(r), neq);
  MapVec jacData(jacVals, nnz);

  rFunc(tres, uu, up, r, userData);
  Eigen::VectorXd u0 = u;
  Eigen::VectorXd r0 = resvec;

  Eigen::VectorXd d(neq), fjacd(neq), fjac(nnz);
  const int col = 1;
  if (alpha) {
    for (int numgrp = 1; numgrp <= maxgrp; numgrp++) {
      for (int j = 0; j < neq; j++) {
        d[j] = 0;
        if (ngrp[j] == numgrp) {
          d[j] = stepLen(u0[j]);
        }
        u[j] = u0[j] + d[j];
      }
      rFunc(tres, uu, up, r, userData);
      fjacd = resvec - r0;
      fdjs_(&neq, &neq, &col, indrow.data(), jpntr.data(), ngrp.data(),
        &numgrp, d.data(), fjacd.data(), fjac.data());
    }
    jacData = alpha * fjac;
  }
  else
    jacData.setZero();

  if (beta) {
    u = u0;
    Eigen::VectorXd up0 = uDot;
    for (int numgrp = 1; numgrp <= maxgrp; numgrp++) {
      for (int j = 0; j < neq; j++) {
        d[j] = 0;
        if (ngrp[j] == numgrp) {
          d[j] = stepLen(up0[j]);
        }
        uDot[j] = up0[j] + d[j];
      }
      rFunc(tres, uu, up, r, userData);
      fjacd = resvec - r0;
      fdjs_(&neq, &neq, &col, indrow.data(), jpntr.data(), ngrp.data(),
        &numgrp, d.data(), fjacd.data(), fjac.data());
    }
    jacData += beta * fjac;
  }

  copyIndices(jacColPtrs, jacRowIndices);
}

void FiniteDiffJacobian::calcJacobianCD(double tres, double alpha,
  double beta, N_Vector uu, N_Vector up, N_Vector r,
  IDAResFn rFunc, void *userData, double *jacVals,
  int *jacColPtrs, int *jacRowIndices)
{
  typedef Eigen::Map<Eigen::VectorXd> MapVec;
  MapVec u(NV_DATA_S(uu), neq);
  MapVec uDot(NV_DATA_S(up), neq);
  MapVec resvec(NV_DATA_S(r), neq);
  MapVec jacData(jacVals, nnz);

  rFunc(tres, uu, up, r, userData);
  Eigen::VectorXd u0 = u;
  Eigen::VectorXd r0 = resvec;

  Eigen::VectorXd d(neq), fjacd(neq), fjac(nnz), rpd(neq);
  SunVector umd(neq);
  const int col = 1;
  if (alpha) {
    for (int numgrp = 1; numgrp <= maxgrp; numgrp++) {
      for (int j = 0; j < neq; j++) {
        d[j] = 0;
        if (ngrp[j] == numgrp) {
          d[j] = stepLen(u0[j]);
        }
        u[j] = u0[j] + d[j];
        umd[j] = u0[j] - d[j];
      }
      rFunc(tres, uu, up, r, userData);
      rpd = resvec;
      rFunc(tres, umd.getNV(), up, r, userData);
      fjacd = (rpd-resvec)/2.;
      fdjs_(&neq, &neq, &col, indrow.data(), jpntr.data(), ngrp.data(),
        &numgrp, d.data(), fjacd.data(), fjac.data());
    }
    jacData = alpha*fjac;
  }
  else
    jacData.setZero();

  if (beta) {
    u = u0;
    Eigen::VectorXd up0 = uDot;
    for (int numgrp = 1; numgrp <= maxgrp; numgrp++) {
      for (int j = 0; j < neq; j++) {
        d[j] = 0;
        if (ngrp[j] == numgrp) {
          d[j] = stepLen(up0[j]);
        }
        uDot[j] = up0[j] + d[j];
        umd[j] = up0[j] - d[j];
      }
      rFunc(tres, uu, up, r, userData);
      rpd = resvec;
      rFunc(tres, uu, umd.getNV(), r, userData);
      fjacd = (rpd - resvec) / 2.;
      fdjs_(&neq, &neq, &col, indrow.data(), jpntr.data(), ngrp.data(),
        &numgrp, d.data(), fjacd.data(), fjac.data());
    }
    jacData += beta*fjac;
  }

  copyIndices(jacColPtrs, jacRowIndices);
}

void FiniteDiffJacobian::copyIndices(int *outerInd, int *innerInd) {
  // copy the row and column pointers
  indrow.array() -= 1;
  jpntr.array() -= 1;
  std::copy_n(jpntr.data(), neq + 1, outerInd);
  std::copy_n(indrow.data(), nnz, innerInd);
  indrow.array() += 1;
  jpntr.array() += 1;
}



