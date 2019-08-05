#include <boost/timer.hpp>

#include "ExampleHeatCond.h"
#include "PDE1dImpl.h"
#include "PDE1dOptions.h"
#include "PDESolution.h"
#include "PDEModel.h"
#include "SunVector.h"

void PDE1dWarningMsg(const char *id, const char *msg) {
  printf("Warning: %s, %s\n", id, msg);
}

int main()
{
  double L=1, tFinal=.05;
  const int nel=11, nt=5;
  boost::timer timer;
  ExampleHeatCond pde(L, nel, tFinal, nt);

#if 1
  ShapeFunctionManager sfm;
  PDEModel model(pde.getMesh(), 1, pde.getNumPDE(), sfm);
  PDESolution pdeSol(pde, model, 1);
  PDE1dOptions opts;
  try {
    opts.setAbsTol(1e-4);
    //opts.setICDiagnostics(3);
    //opts.setJacDiagnostics(1);
    PDE1dImpl pdeImpl(pde, opts);
    int err = pdeImpl.solveTransient(pdeSol);
    if (err)
      return 1;
  }
  catch (const std::exception &ex) {
    printf("exception caught: %s\n", ex.what());
  }
  catch (...) {
    printf("unknown exception caught\n");
  }
  const RealMatrix &u = pdeSol.getSolution();
  printf("size u=%d,%d\n", u.rows(), u.cols());
  auto ntr = u.rows();
  FILE *fp = fopen("pde1D.out", "w");
  RealVector x = pde.getMesh();
  auto nx = x.size();
  int neq = pde.getNumPDE();
  int ii = 0;
  for (int i = 0; i < nx; i++) {
    fprintf(fp, " %15.6e", x(i));
    for (int j = 0; j < neq; j++)
      fprintf(fp, " %15.6e", u(ntr-1, ii+j));
    fprintf(fp, "\n");
    ii += neq;
  }
  fclose(fp);

#elif 1
  // plot IC
  RealVector x = pde.getMesh();
  int nx = x.size();
  int neq = pde.getNumPDE();
  RealVector ic(neq);
  FILE *fp = fopen("pde1D.out", "w");
  for (int i = 0; i < nx; i++) {
    pde.evalIC(x(i), ic);
    fprintf(fp, " %15.6e", x(i));
    for (int j = 0; j < neq; j++)
      fprintf(fp, " %15.6e", ic(j));
    fprintf(fp, "\n");
  }
  fclose(fp);
#elif 0
  PDE1dOptions opts;
  PDE1dImpl pdeImpl(pde, opts);
  pdeImpl.testMats();
#endif
  printf("Elapsed time = %7.3f seconds.\n", timer.elapsed());
}