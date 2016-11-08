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

#include <unsupported/Eigen/CXX11/Tensor>
typedef Eigen::Tensor<double, 3> RealTensor3;


#endif /* PDE1DLIB_MATRIXTYPES_H_ */
