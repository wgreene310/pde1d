#pragma once

#include "PDE1dException.h"

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_sparse.h>  
#include <sundials/sundials_linearsolver.h>
#include <nvector/nvector_serial.h> 

#include <Eigen/SparseCore>
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, sunindextype> SSparseMat;
typedef Eigen::Map<SSparseMat> ESM;

class SunSparseMapX : public ESM {
public:
  SunSparseMapX(SUNMatrix a) :
    ESM(SUNSparseMatrix_Rows(a), SUNSparseMatrix_Columns(a), SUNSparseMatrix_NNZ(a),
      SUNSparseMatrix_IndexPointers(a), SUNSparseMatrix_IndexValues(a),
      SUNSparseMatrix_Data(a)) {

  }
};

#include<Eigen/SparseLU>

class EigenSUNSparseSolver : public _generic_SUNLinearSolver {
public:
  typedef EigenSUNSparseSolver EigenLU;
  EigenSUNSparseSolver(N_Vector y, SUNMatrix A) {
    lastFlag = 0;
    ops = &_ops;
    /* Attach operations */
    ops->gettype = LinSolGetType;
    ops->initialize = LinSolInitialize;
    ops->setup = LinSolSetup;
    ops->solve = LinSolSolve;
    ops->lastflag = LinSolLastFlag;
    ops->space = LinSolSpace;
    ops->free = LinSolFree;
    ops->setatimes = NULL;
    ops->setpreconditioner = NULL;
    ops->setscalingvectors = NULL;
    ops->numiters = NULL;
    ops->resnorm = NULL;
    ops->resid = NULL;
  }
private:
  inline static EigenLU &getSolver(SUNLinearSolver S) {
    return *static_cast<EigenSUNSparseSolver*>(S);
  }
  static SUNLinearSolver_Type LinSolGetType(SUNLinearSolver S) {
    return(SUNLINEARSOLVER_DIRECT);
  }
  static int LinSolInitialize(SUNLinearSolver S) {
    return 0;
  }
  static int LinSolSetup(SUNLinearSolver S, SUNMatrix A) {
    auto &s = getSolver(S);
    SunSparseMapX a(A);
    //cout << "nnz=" << a.nonZeros() << endl;
    //cout << a << endl;
    //a.prune(1.0); // remove zero entries for eigen
    s.solver.compute(a);
    s.lastFlag = s.solver.info();
    if (s.lastFlag != Eigen::Success) {
      // decomposition failed
      return 1;
    }
    return 0;
  }
  static int LinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
    N_Vector b, realtype tol) {
    auto &s = getSolver(S);
    Eigen::Map<Eigen::VectorXd> xx(NV_DATA_S(x), NV_LENGTH_S(x));
    Eigen::Map<Eigen::VectorXd> bb(NV_DATA_S(b), NV_LENGTH_S(b));
    xx = s.solver.solve(bb);
    if (s.solver.info() != Eigen::Success) {
      // solving failed
      return 1;
    }
    return 0;
  }
  static long int LinSolLastFlag(SUNLinearSolver S) {
    auto &s = getSolver(S);
    return s.lastFlag;
  }
  static int LinSolSpace(SUNLinearSolver S, long int *lenrwLS,
    long int *leniwLS) {
    auto &s = getSolver(S);
    throw PDE1dException("pde1d:solver_internal", "LinSolSpace called.");
    return 0;
  }
  static int LinSolFree(SUNLinearSolver S) {
    auto &s = getSolver(S);
    return 0;
  }
  _generic_SUNLinearSolver_Ops _ops;
  Eigen::SparseLU<ESM> solver;
  int lastFlag;
};