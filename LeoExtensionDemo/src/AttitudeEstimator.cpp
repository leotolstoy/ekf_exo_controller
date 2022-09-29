/*
 * AttitudeEstimator.cpp
 *
 *    Created on: May 13, 2014
 *        Author: Gray Thomas
 *   Affiliation: Human Centered Robotics Lab, University of Texas at Austin
 *
 */

#include "AttitudeEstimator.hpp"
#include <SpatialTransform.hpp>
#include <iostream>

namespace gu
{
namespace attitudeEstimator
{

AttitudeEstimator::AttitudeEstimator() :
			center(),
			directionCosines(),
			rawDirectionCosines(),
//			directionQuaternion(),
			directionQuat4(),
			directionTest(),
			regressor(),
			symmetricRegressionProduct(),
			symmetricRegressionProductInverse(),
			linearSolutionMatrix(),
			observations(),
			constellation(),
			numPoints(7)
{
	center.setZero(3);
	directionCosines.setZero(3, 3);
//	directionQuaternion.setZero(4);
	directionQuat4.setIdentity();
	directionTest.setZero(4, 4);
	symmetricRegressionProduct.setZero(12, 12);
	symmetricRegressionProductInverse.setZero(12, 12);
	linearSolutionMatrix.setZero(12, 3 * numPoints);
	regressor.setZero(3 * numPoints, 12);
	observations.setZero(3 * numPoints);
	constellation.setZero(3, numPoints);
	constellation << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
	populateVectorXdFromMatrix3Xd(observations, constellation);
//	std::cout << constellation << std::endl;
//	std::cout << observations << std::endl;
}

void AttitudeEstimator::solveDirectly(Matrix3Xd observedPoints)
{
	populateVectorXdFromMatrix3Xd(observations, observedPoints);
	Vector12d result = linearSolutionMatrix * observations;
	center << result.segment(0, 3);
	rawDirectionCosines << result.segment(3, 3).transpose(), result.segment(6, 3).transpose(), result.segment(9, 3).transpose();
	gu::vectorGeometry::closestQuaternionToMatrix(directionQuat4, rawDirectionCosines);
	gu::vectorGeometry::directionCosinesFromQuaternion(directionQuat4, directionCosines);
}
void moveConstellationCentroidToOrigin(gu::attitudeEstimator::Matrix3Xd& constellation)
{
	Vector3d originShift = constellation.rowwise().sum() * (1.0 / 7.0);
	constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);
}
Vector3d svdConstellation(Matrix3Xd& constellation)
{
	Matrix3d constellationMatrix;
	constellationMatrix.setZero();
	for (int i = 0; i < constellation.cols(); i++)
		constellationMatrix += constellation.col(i) * constellation.col(i).transpose();
	JacobiSVD<Matrix3d> svdy(constellationMatrix);
	return Vector3d(svdy.singularValues());
}
void AttitudeEstimator::setConstellation(Matrix3Xd constellation)
{
	numPoints = constellation.cols();
	regressor.setZero(3 * numPoints, 12);
	observations.setZero(3 * numPoints);
	linearSolutionMatrix.setZero(12, 3 * numPoints);
	this->constellation.resize(3, numPoints);
	this->constellation << constellation;
	std::cout << "hey" << std::endl;
	for (int i = 0; i < numPoints; i++)
	{
		regressor.block(3 * i, 0, 3, 3) << Matrix3d::Identity();
		regressor.block(3 * i, 3, 3, 3) << constellation.col(i).transpose(), MatrixXd::Zero(2, 3);
		regressor.block(3 * i, 6, 3, 3) << MatrixXd::Zero(1, 3), constellation.col(i).transpose(), MatrixXd::Zero(1,
				3);
		regressor.block(3 * i, 9, 3, 3) << MatrixXd::Zero(2, 3), constellation.col(i).transpose();
	}
	this->symmetricRegressionProduct = regressor.transpose() * regressor;
	this->symmetricRegressionProductInverse = symmetricRegressionProduct.inverse();
	linearSolutionMatrix << symmetricRegressionProductInverse * regressor.transpose();

//	std::cout << "svd of constellation for rotation:\n" << svdConstellation(constellation) << std::endl;
}

void AttitudeEstimator::setConstellationWithPrior(const Matrix3Xd& constellation, double priorPositionWeight,
													double priorOrientationWeight)
{
	numPoints = constellation.cols();
	regressor.setZero(3 * numPoints + 12, 12);
	observations.setZero(3 * numPoints + 12);
	linearSolutionMatrix.setZero(12, 3 * numPoints + 12);
	this->constellation.resize(3, numPoints);
	this->constellation << constellation;
	for (int i = 0; i < numPoints; i++)
	{
		regressor.block(3 * i, 0, 3, 3) << Matrix3d::Identity();
		regressor.block(3 * i, 3, 3, 3) << constellation.col(i).transpose(), MatrixXd::Zero(2, 3);
		regressor.block(3 * i, 6, 3, 3) << MatrixXd::Zero(1, 3), constellation.col(i).transpose(), MatrixXd::Zero(1,
				3);
		regressor.block(3 * i, 9, 3, 3) << MatrixXd::Zero(2, 3), constellation.col(i).transpose();
	}
	regressor.block(numPoints * 3, 0, 3, 3) << Matrix3d::Identity();
	regressor.block(numPoints * 3 + 3, 3, 9, 9) << MatrixXd::Identity(9, 9);
	VectorXd weights(numPoints * 3 + 12);
	weights.setOnes();
	weights.tail(9) *= priorOrientationWeight;
	weights.segment(3 * numPoints, 3) *= priorPositionWeight;
	this->symmetricRegressionProduct = regressor.transpose() * weights.asDiagonal() * regressor;
	this->symmetricRegressionProductInverse = symmetricRegressionProduct.inverse();
	linearSolutionMatrix << symmetricRegressionProductInverse * regressor.transpose() * weights.asDiagonal();
//	std::cout << "regressor with priors\n" << regressor << std::endl;
//	std::cout << "solutionMatrix with priors\n" << linearSolutionMatrix << std::endl;
}

void AttitudeEstimator::solveDirectlyWithPrior(const Matrix3Xd& observedPoints, const Vector3d& priorPosition,
												const Matrix3d& priorOrientation)
{
	unsigned numPoints = observedPoints.cols();

	observations.resize(3*numPoints+12);
	populateVectorXdFromMatrix3Xd(observations, observedPoints);
	observations.segment(3 * numPoints, 3) << priorPosition;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			observations(3 * numPoints + 3 + i * 3 + j) = priorOrientation(i, j);
	Vector12d result = linearSolutionMatrix * observations;
	center << result.segment(0, 3);
	rawDirectionCosines << result.segment(3, 3).transpose(), result.segment(6, 3).transpose(), result.segment(9, 3).transpose();
	gu::vectorGeometry::closestQuaternionToMatrix(directionQuat4, rawDirectionCosines);
	gu::vectorGeometry::directionCosinesFromQuaternion(directionQuat4, directionCosines);
}

Quaterniond AttitudeEstimator::getQuaternion()
{
	Quaterniond ret(directionQuat4);
	return ret;
}
Vector3d AttitudeEstimator::getPoint()
{
	Vector3d ret(center);
	return ret;
}
Matrix3d AttitudeEstimator::getDCM()
{
	Matrix3d ret(directionCosines);
	return ret;
}
void sayHi(void)
{
	std::cout << "Hello!" << std::endl;
}
void populateVectorXdFromMatrix3Xd(VectorXd& vec, const Matrix3Xd& mat)
{
	for (int i = 0; i < mat.cols(); i++)
		vec.segment(3 * i, 3) << mat.col(i);
}

}


}
