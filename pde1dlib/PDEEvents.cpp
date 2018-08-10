#include "PDEEvents.h"
#include <PDE1dDefn.h>
#include <PDEModel.h>



PDEEvents::PDEEvents(PDE1dDefn &pde, PDEModel &pdeModel) :
  pde(pde), pdeModel(pdeModel)
{
  numEvents = pde.getNumEvents();
  if (numEvents) {
    size_t nn = pde.getMesh().size();
    int numPde = pde.getNumPDE();
    eventsVal.resize(numEvents);
    eventsIsTerminal.resize(numEvents);
    eventsDirection.resize(numEvents);
    eventsDirection.setZero();
    eventsU.resize(numPde, nn);
  }
}

void PDEEvents::calcEvents(double time, const RealVector &u, double *gOut)
{
  MapVec g(gOut, numEvents);
  pdeModel.globalToMeshVec(u, eventsU);
  pde.evalEvents(time, eventsU, eventsVal, eventsIsTerminal,
    eventsDirection);
  g = eventsVal;
}

bool PDEEvents::isTerminalEvent(double time, const RealVector &u,
  const IntVector &eventsFound) {
  bool doTerm = false;
  pdeModel.globalToMeshVec(u, eventsU);
  pde.evalEvents(time, eventsU, eventsVal, eventsIsTerminal,
    eventsDirection);
  for (int j = 0; j < numEvents; j++) {
    int jFlag = eventsFound(j);
#if 0
    printf("Event %d: jFlag=%d, eventsVal=%12.3e, "
      "eventsIsTerminal=%g, eventsDirection=%g\n",
      j, jFlag, eventsVal(j), eventsIsTerminal(j),
      eventsDirection(j));
#endif
    if (jFlag && eventsIsTerminal(j)) {
      //printf("terminating on event=%d, tret=%12.4e\n", j, time);
      doTerm = true;
    }
  }
  return doTerm;
}

PDEEvents::~PDEEvents()
{
}
