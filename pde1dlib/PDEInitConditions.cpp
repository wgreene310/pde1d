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

/*
 * PDEInitConditions.cpp
 *
 *  Created on: Sep 27, 2016
 *      Author: bgreene
 */

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
#if 0
#define EIGEN_USE_LAPACKE 1
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif

#include <Eigen/QR>

#include <ida/ida.h>

#include "PDEInitConditions.h"
#include "PDE1dImpl.h"
#include "SunVector.h"
#include "PDE1dOptions.h"
#include "PDE1dException.h"
#include "PDE1dWarningMsg.h"
#include "util.h"

#define USE_LAPACK_QR 0

#if USE_LAPACK_QR
extern "C" {
  void dgeqp3_(const int *M, const int *N, double *A,
    const int *LDA, int *JPVT, double *TAU, double *WORK,
    const int *LWORK, int *INFO);
  void dormqr_(const char *SIDE, const char *TRANS, const int *M, const int *N, const int *K, double *A,
    const int *LDA, double *TAU, double *C, const int *LDC, double *WORK,
    int *LWORK, int &INFO);
  void dtrsm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG, const int *M, const int *N,
    const double *ALPHA, double *A, const int *LDA, double *B, const int *LDB);
}
#endif


PDEInitConditions::PDEInitConditions(void *idaMem, PDE1dImpl &pdeImpl,
  const SunVector &u0, const SunVector &up0) :
  idaMem(idaMem), pdeImpl(pdeImpl), u0(u0), up0(up0)
{
  icMeth = pdeImpl.getOptions().getICMethod();
  size_t ne = u0.rows();
  u0C = std::unique_ptr<SunVector>(new SunVector(ne));
  up0C = std::unique_ptr<SunVector>(new SunVector(ne));

}

PDEInitConditions::ICPair PDEInitConditions::init() {
  ICPair icPair;
  const RealVector &tspan = pdeImpl.getTimePoints();
  double t0 = tspan[0];
  if (!icMeth) {
    calcShampineAlgo(t0, *u0C, *up0C);
    icPair.first = u0C.get();
    icPair.second = up0C.get();
  }
  else {
    SunVector &u0x = const_cast<SunVector&>(u0);
    SunVector &up0x = const_cast<SunVector&>(up0);
    icPair.first = &u0x;
    icPair.second = &up0x;
  }
  return icPair;
}

void PDEInitConditions::calcShampineAlgo(double t0,
  SunVector &yNew, SunVector &ypNew)
{
  const size_t numEqns = u0.rows();
  SparseMat dfDy(numEqns, numEqns), dfDyp(numEqns, numEqns);
  Eigen::MatrixXd dfDypDense(numEqns, numEqns);
  const int maxIter = 10;
  RealVector dy(numEqns), dyp(numEqns), dypP(numEqns);
  SunVector res(numEqns);
  yNew = u0;
  ypNew = up0;
  int it = 0;
  bool converged = false;
  double maxRes = 1;
  int diag = pdeImpl.getOptions().getICDiagnostics();
  if (diag > 1) {
    cout << "u0=" << yNew.transpose() << endl;
    cout << "up0=" << ypNew.transpose() << endl;
  }
  double absTol = pdeImpl.getOptions().getAbsTol();
  Eigen::ColPivHouseholderQR<RealMatrix> qr;
  while (it++ < maxIter) {
    pdeImpl.calcRHSODE(t0, yNew, ypNew, res);
#if 1
    // rms tolerance
    double resRms = sqrt(res.dot(res)) / (double)numEqns;
    if (diag)
      printf("IC: iter = %d, resRms=%12.3e\n", it, resRms);
    if (resRms < absTol) {
      converged = true;
      break;
    }
#else
    // max tolerance
    auto absRes = res.array().abs();
    maxRes = absRes.maxCoeff();
    printf("iter = %d, maxRes=%12.3e\n", it, maxRes);
    if ((absRes < options.getAbsTol()).all()) {
      converged = true;
      break;
    }
#endif
    if (diag > 1)
      cout << "res=" << res.transpose() << endl;
    pdeImpl.calcJacobian(t0, 1, 0, yNew, ypNew, res, dfDy);
    pdeImpl.calcJacobian(t0, 0, 1, yNew, ypNew, res, dfDyp);
    if (diag > 2) {
      cout << "dfDy\n" << dfDy.toDense() << endl;
      cout << "dfDyp\n" << dfDyp.toDense() << endl;
    }
#if 0
    if (it == 1) {
      Eigen::saveMarket(dfDyp, "dfDyp.mtx");
    }
#endif
#if 0
    SparseQR qr(dfDyp);
#else
    dfDypDense = dfDyp.toDense();
    qr.compute(dfDypDense);
#endif

    RealVector d = -(qr.matrixQ().transpose()*res);
    RealMatrix S = qr.matrixQ().transpose()*dfDy.toDense();
    size_t rnk = qr.rank();
    size_t numAlgVars = numEqns - rnk;
    if(diag)
      printf("IC: numEqns=%zd, rnk=%zd, numAlgVars=%zd\n",
      numEqns, rnk, numAlgVars);

#if 0
    Eigen::VectorXd rDiag = qr.matrixR().diagonal();
    cout << "Rdiag=" << rDiag.transpose() << endl;
#endif

    auto S2122 = S.bottomRows(numAlgVars);
    //cout << S2122 << endl;
    Eigen::ColPivHouseholderQR<RealMatrix> qrS(S2122);
    if (diag)
      cout << "IC: rank of S=" << qrS.rank() << endl;
    if (qrS.rank() != numAlgVars) {
      cout << "Error detected in computation of consistent initial conditions." << endl;
    }
    dy = qrS.solve(d.bottomRows(numAlgVars));
    if (diag>1)
      cout << "dy=" << dy.transpose() << endl;
    auto R11 = qr.matrixR().topLeftCorner(rnk, rnk).triangularView<Eigen::Upper>();
    RealVector w1p = R11.solve(d.topRows(rnk) - S.topRows(rnk)*dy);
    dypP.setZero();
    dypP.topRows(rnk) = w1p;
    dyp = qr.colsPermutation()*dypP;
    if (diag>1)
      cout << "dyp=" << dyp.transpose() << endl;

#if USE_LAPACK_QR
    const int nInt = (int)numEqns;
    int info = 0;
    Eigen::VectorXi jpiv(nInt);
    jpiv.setZero();
    const int nb = 64;
    int lwork = static_cast<int>(2 * numEqns + (numEqns + 1)*nb);
    Eigen::VectorXd tau(numEqns), work(lwork);
#if 1
    dgeqp3_(&nInt, &nInt, dfDypDense.data(), &nInt, jpiv.data(),
      tau.data(), work.data(), &lwork, &info);
#endif
    cout << "info=" << info << endl;
    cout << "jpiv=" << jpiv.transpose() << endl;
    Eigen::VectorXd rii = dfDypDense.diagonal();
    cout << "rii=" << rii.transpose() << endl;
    double tol = absTol*abs(rii(0));
    int rankLap = 1;
    for (int i = 1; i < numEqns; i++) {
      if (abs(rii(i)) < tol) break;
      rankLap++;
    }
    cout << "rankLap=" << rankLap << endl;
    int nrhs = 1;
#if 0
    dormqr_("L", "T", &nInt, &nrhs, &nInt, dfDypDense.data(), &nInt, tau.data(),
      b.data(), &m, work.data(), &lwork, info);
#endif
    cout << "info=" << info << endl;
#endif

    yNew += dy;
    ypNew += dyp;
  }

  if (!converged) {
    printf("Unable to obtain a consistent set of initial conditions.\n"
      "Maximum error in the residual is %12.3e.\n", maxRes);
  }
}

void PDEInitConditions::calcSundialsAlgo(double tf,
  SunVector &yNew, SunVector &ypNew)
{
#if DAE_Y_INIT
  int ier = IDACalcIC(idaMem, IDA_Y_INIT, tf);
#else
  int ier = IDACalcIC(idaMem, IDA_YA_YDP_INIT, tf);
#endif
  if (ier < 0)
    icFailErr();
  ier = IDAGetConsistentIC(idaMem, yNew.getNV(), ypNew.getNV());
  if (ier < 0)
    throw PDE1dException("pde1d:get_consistent_ic",
    "Unable to retrieve consistent initial conditions from IDA.");
}

void PDEInitConditions::icFailErr()
{
  throw PDE1dException("pde1d:consistent_ic",
    "Unable to calculate a consistent initial solution to the PDE.\n"
    "Often this is caused by an incorrect specification of the boundary conditions.");
}

void PDEInitConditions::update()
{
  const RealVector &tspan = pdeImpl.getTimePoints();
  if (icMeth > 0) {
    double tf = tspan.tail<1>()[0];
    int ier;
    if (icMeth == 1) {
      ier = IDACalcIC(idaMem, IDA_YA_YDP_INIT, tf);
    }
    else if (icMeth == 2) {
      ier = IDACalcIC(idaMem, IDA_Y_INIT, tf);
    }
    else
      throw PDE1dException("pde1d:consistent_ic_invalid_method",
      "Invalid method for calculation of consistent initial conditions.");
    if (ier < 0)
      icFailErr();
    ier = IDAGetConsistentIC(idaMem, u0C->getNV(), up0C->getNV());
    if (ier < 0)
      throw PDE1dException("pde1d:get_consistent_ic",
      "Unable to retrieve consistent initial conditions from IDA.");
  }

  int diag = pdeImpl.getOptions().getICDiagnostics();
  if (diag) {
    //cout << "IC: yNew=" << yNew.transpose() << endl;
    //cout << "IC: ypNew=" << ypNew.transpose() << endl;
    ::print(u0C->getNV(), "yNew");
    ::print(up0C->getNV(), "ypNew");
  }
}

void PDEInitConditions::compareInitConditions()
{
  double atol = pdeImpl.getOptions().getAbsTol();
  double maxDiff = (u0 - *u0C).cwiseAbs().maxCoeff();
  if (maxDiff > 100 * atol)
    PDE1dWarningMsg("pde1d:init_cond_changed",
    "User-defined initial conditions were changed "
    "to create a consistent solution to the equations at "
    "the initial time.\n");
}

void PDEInitConditions::print() const
{
  cout << "Initial y=" << u0C->transpose() << endl;
  cout << "Initial y_dot=" << up0C->transpose() << endl;
}

