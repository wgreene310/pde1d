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

#include <string>
#include <vector>

#include <Eigen/Core>

#include <mex.h>

class MexInterface {
public:
  typedef Eigen::VectorXd RealVector;
  MexInterface(int maxOutArgs=1);
  ~MexInterface();
  void print(const mxArray *a, const char *name);
  std::string getFuncNameFromHandle(const mxArray *fh);
  void callMatlab(const mxArray *inArgs[], int nargin,
    RealVector *outArgs[], int nargout);
  template<class T>
  static mxArray *toMxArray(const T &a)
  {
    mxArray *ma = mxCreateDoubleMatrix(a.rows(), a.cols(), mxREAL);
    double *dest = mxGetPr(ma);
    const typename T::Scalar *src = a.data();
    std::copy_n(src, a.cols()*a.rows(), dest);
    return ma;
  }
  static Eigen::MatrixXd fromMxArray(const mxArray *a);
  static Eigen::VectorXd fromMxArrayVec(const mxArray *a);
private:
  const int maxOutArgs;
  std::vector<mxArray*> matOutArgs;
};

