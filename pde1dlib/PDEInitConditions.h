/*
 * PDEInitConditions.h
 *
 *  Created on: Sep 27, 2016
 *      Author: bgreene
 */

#ifndef PDE1DLIB_PDEINITCONDITIONS_H_
#define PDE1DLIB_PDEINITCONDITIONS_H_

#include <memory>

class SunVector;
class PDE1dImpl;

class PDEInitConditions {
public:
  PDEInitConditions(void *idaMem, PDE1dImpl &pdeImpl,
    const SunVector &u0, const SunVector &up0);
  typedef std::pair<SunVector*, SunVector*> ICPair;
  ICPair init();
  void update();
  const SunVector &getU0() {
    return *u0C.get();
  }
private:
  void calcShampineAlgo(double t0,
    SunVector &yNew, SunVector &ypNew);
  void calcSundialsAlgo(double tf,
    SunVector &yNew, SunVector &ypNew);
  void icFailErr();
  bool compareInitConditions();
  int icMeth;
  void *idaMem;
  PDE1dImpl &pdeImpl;
  const SunVector &u0, &up0;
  std::unique_ptr<SunVector> u0C, up0C;
};

#endif /* PDE1DLIB_PDEINITCONDITIONS_H_ */
