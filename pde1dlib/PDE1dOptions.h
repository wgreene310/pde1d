#pragma once
class PDE1dOptions
{
public:
  PDE1dOptions(double relTol = 1e-3, double absTol = 1e-6);
  double getRelTol() const { return relTol;  }
  double getAbsTol() const { return absTol;  }
  void setRelTol(double tol) { relTol = tol; }
  void setAbsTol(double tol) { absTol = tol; }
private:
  double relTol, absTol;
};

