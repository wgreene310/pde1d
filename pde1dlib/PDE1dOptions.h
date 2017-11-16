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

#ifndef PDE1dOptions_h
#define PDE1dOptions_h

class PDE1dOptions
{
public:
  PDE1dOptions(double relTol = 1e-3, double absTol = 1e-6) :
    relTol(relTol), absTol(absTol) {
    vectorizedFuncs = stats = false;
    maxSteps = 10000;
    ICMethod = 0;
    ICDiagnostics = 0;
    polyOrder = 1;
    viewMesh = 1;
    useDiagMassMat = false;
  }
  double getRelTol() const { return relTol;  }
  double getAbsTol() const { return absTol;  }
  void setRelTol(double tol) { relTol = tol; }
  void setAbsTol(double tol) { absTol = tol; }
  bool isVectorized() const { return vectorizedFuncs; }
  void setVectorized(bool isVec) { vectorizedFuncs = isVec; }
  int getMaxSteps() const { return maxSteps;  }
  void setMaxSteps(int maxSteps) { this->maxSteps = maxSteps;  }
  bool printStats() const { return stats; }
  void setPrintStats(bool sts) { stats = sts; }
  int getICMethod() const { return ICMethod; }
  void setICMethod(int meth) { ICMethod = meth;  }
  int getICDiagnostics() const { return ICDiagnostics;  }
  void setICDiagnostics(int diag) { ICDiagnostics = diag; }
  void setPolyOrder(int pOrder) { polyOrder = pOrder;  }
  int getPolyOrder() const {
    return polyOrder;
  } 
  void setViewMesh(int intVUMesh) { viewMesh = intVUMesh; }
  int getViewMesh() const {
    return viewMesh;
  }
  void setDiagMassMat(bool useDiagMassMat) {
    this->useDiagMassMat = useDiagMassMat;
  }
  bool getDiagMassMat() const {
    return useDiagMassMat;
  }
private:
  double relTol, absTol;
  bool vectorizedFuncs;
  int maxSteps;
  bool stats;
  int ICMethod;
  int ICDiagnostics;
  int polyOrder;
  int viewMesh;
  bool useDiagMassMat;
};

#endif

