#pragma once

#include <nvector/nvector_serial.h>
#include <ida/ida.h>
#include <ida/ida_klu.h>

#include <Eigen/SparseCore>

struct SunSparseMat : public _SlsMat {
  SunSparseMat(Eigen::SparseMatrix<double> &eMat) {
    M = eMat.rows();
    N = eMat.cols();
    NNZ = eMat.nonZeros();
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
#if 0
  void calcJacobian(SparseMat &jac, IDAResFn rf, double tres, double alpha,
    void *userData);
  void calcJacobian(SlsMat jac, IDAResFn rf, double tres, double alpha,
    void *userData);
#else
  void calcJacobian(double tres, double alpha,
    N_Vector uu, N_Vector up, N_Vector r,
    IDAResFn rf, void *userData, SparseMat &Jac);
  void calcJacobian(double tres, double alpha,
    N_Vector uu, N_Vector up, N_Vector r,
    IDAResFn rf, void *userData, SlsMat Jac);
#endif
private:
  int neq, nnz;
  Eigen::VectorXi indrow, jpntr, ngrp;
  int maxgrp, mingrp;
  double sqrtEps;
};

