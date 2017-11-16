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
  const SunVector &getUp0() {
    return *up0C.get();
  }
  void print() const;
private:
  void calcShampineAlgo(double t0,
    SunVector &yNew, SunVector &ypNew);
  void calcSundialsAlgo(double tf,
    SunVector &yNew, SunVector &ypNew);
  void icFailErr();
  void compareInitConditions();
  int icMeth;
  void *idaMem;
  PDE1dImpl &pdeImpl;
  const SunVector &u0, &up0;
  std::unique_ptr<SunVector> u0C, up0C;
};

#endif /* PDE1DLIB_PDEINITCONDITIONS_H_ */
