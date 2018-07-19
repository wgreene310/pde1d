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
 * MatrixTypes.h
 *
 *  Created on: Sep 5, 2016
 *      Author: bgreene
 */

#ifndef PDE1DLIB_MATRIXTYPES_H_
#define PDE1DLIB_MATRIXTYPES_H_

#include <Eigen/Core>

typedef Eigen::VectorXd RealVector;
typedef Eigen::MatrixXd RealMatrix;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::VectorXi IntVector;

#include <unsupported/Eigen/CXX11/Tensor>
typedef Eigen::Tensor<double, 3> RealTensor3;


#endif /* PDE1DLIB_MATRIXTYPES_H_ */
