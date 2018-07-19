#pragma once

#include "MatrixTypes.h"

class PDE1dDefn;
class PDEModel;

class PDEEvents
{
public:
  PDEEvents(PDE1dDefn &pde, PDEModel &pdeModel);
  void calcEvents(double time, const RealVector &u, double *gOut);
  bool isTerminalEvent(double time, const RealVector &u,
    const IntVector &eventsFound);
  ~PDEEvents();
private:
  int numEvents;
  RealVector eventsVal, eventsIsTerminal, eventsDirection;
  RealMatrix eventsU;
  PDE1dDefn &pde;
  const PDEModel &pdeModel;
};
