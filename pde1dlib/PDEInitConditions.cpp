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

#include <Eigen/QR>

#include <ida/ida.h>

#include "PDEInitConditions.h"
#include "PDE1dImpl.h"
#include "SunVector.h"
#include "PDE1dOptions.h"
#include "PDE1dException.h"
#include "PDE1dWarningMsg.h"
#include "util.h"


PDEInitConditions::PDEInitConditions(void *idaMem, PDE1dImpl &pdeImpl,
  const SunVector &u0, const SunVector &up0) :
  idaMem(idaMem), pdeImpl(pdeImpl), u0(u0), up0(up0)
{
  icMeth = pdeImpl.getOptions().getICMethod();
  int ne = u0.rows();
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
  const int numEqns = u0.rows();
  SparseMat dfDy(numEqns, numEqns), dfDyp(numEqns, numEqns);
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
  while (it++ < maxIter) {
    pdeImpl.calcRHSODE(t0, yNew, ypNew, res);
#if 1
    // rms tolerance
    double resRms = sqrt(res.dot(res)) / (double)numEqns;
    if (diag)
      printf("IC: iter = %d, resRms=%12.3e\n", it, resRms);
    if (resRms < pdeImpl.getOptions().getAbsTol()) {
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
    //cout << "dfDy\n" << dfDy.toDense() << endl;
    pdeImpl.calcJacobian(t0, 0, 1, yNew, ypNew, res, dfDyp);
    //cout << "dfDyp\n" << dfDyp.toDense() << endl;
#if 0
    if (it == 1) {
      Eigen::saveMarket(dfDyp, "dfDyp.mtx");
    }
#endif
#if 0
    SparseQR qr(dfDyp);
#else
    Eigen::ColPivHouseholderQR<RealMatrix> qr(dfDyp.toDense());
#endif
    RealVector d = -(qr.matrixQ().transpose()*res);
    RealMatrix S = qr.matrixQ().transpose()*dfDy.toDense();
    int rnk = qr.rank();
    int numAlgVars = numEqns - rnk;
    if(diag)
      printf("IC: numEqns=%d, rnk=%d, numAlgVars=%d\n",
      numEqns, rnk, numAlgVars);

#if 0
    Eigen::VectorXd rDiag = qr.matrixR().diagonal();
    cout << "Rdiag=" << rDiag.transpose() << endl;
#endif

    auto S2122 = S.bottomRows(numAlgVars);
    //cout << S2122 << endl;
    Eigen::ColPivHouseholderQR<RealMatrix> qrS(S2122);
    if (qrS.rank() != numAlgVars) {
      printf("Error detected in computation of consistent initial conditions.\n");
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

