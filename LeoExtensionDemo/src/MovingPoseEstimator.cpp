/*
 * AttitudeEstimator.cpp
 *
 *    Created on: May 13, 2014
 *        Author: Gray Thomas
 *   Affiliation: Human Centered Robotics Lab, University of Texas at Austin
 *
 */

#include "MovingPoseEstimator.hpp"
#include "AttitudeEstimator.hpp"
#include <SpatialTransform.hpp>
#include <iostream>

namespace gu
{
namespace attitudeEstimator
{

MovingPoseEstimator::MovingPoseEstimator() :
			center(),
			velocityDeviation(),
			directionCosines(),
			rawDirectionCosines(),
			directionQuaternion(),
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
	directionQuaternion.setZero(4);
	directionTest.setZero(4, 4);
	symmetricRegressionProduct.setZero(12, 12);
	symmetricRegressionProductInverse.setZero(12, 12);
	linearSolutionMatrix.setZero(12, 3 * numPoints);
	regressor.setZero(3 * numPoints, 12);
	observations.setZero(3 * numPoints);
	constellation.setZero(3, numPoints);
	constellation << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
	gu::attitudeEstimator::sayHi();
	populateVectorXdFromMatrix3Xd(observations, constellation);
//	std::cout << constellation << std::endl;
//	std::cout << observations << std::endl;
}

//void MovingPoseEstimator::solveDirectly(Matrix3Xd observedPoints)
//{
//	populateVectorXdFromMatrix3Xd(observations, observedPoints);
//	Vector12d result = linearSolutionMatrix * observations;
//	center << result.segment(0, 3);
//	rawDirectionCosines << result.segment(3, 3).transpose(), result.segment(6, 3).transpose(), result.segment(9, 3).transpose();
////	std::cout << "rawDirectionCosines \n" << rawDirectionCosines << "\ncenter\n" << center << std::endl;
//	Matrix3d& d = rawDirectionCosines;
//	double dxx = d(0, 0), dxy = d(0, 1), dxz = d(0, 2), dyx = d(1, 0), dyy = d(1, 1), dyz = d(1, 2), dzx = d(2, 0),
//			dzy = d(2, 1), dzz = d(2, 2);
//	directionTest.row(0) << (dxx - dyy - dzz), (dyx + dxy), (dzx + dxz), (dyz - dzy);
//	directionTest.row(1) << (dyx + dxy), (dyy - dxx - dzz), (dzy + dyz), (dzx - dxz);
//	directionTest.row(2) << (dzx + dxz), (dzy + dyz), (dzz - dxx - dyy), (dxy - dyx);
//	directionTest.row(3) << (dyz - dzy), (dzx - dxz), (dxy - dyx), (dxx + dyy + dzz);
//	directionTest *= (1.0 / 3.0);
////	std::cout << "directionTest=\n" << directionTest << std::endl;
//	SelfAdjointEigenSolver<Matrix4d> solver(directionTest);
////	std::cout << "eigenvalue\n" << solver.eigenvalues()(3) << std::endl;
//	directionQuaternion << solver.eigenvectors().col(3);
//	directionQuaternion.segment(0, 3) *= -1; // flip direction of transform (conjugate quaternion)
//	if (directionQuaternion(3) < 0)
//		directionQuaternion *= -1;
////	std::cout << "directionQuaternion\n" << directionQuaternion << std::endl;
//	directionQuat4.x() = directionQuaternion(0);
//	directionQuat4.y() = directionQuaternion(1);
//	directionQuat4.z() = directionQuaternion(2);
//	directionQuat4.w() = directionQuaternion(3);
//	directionCosinesFromQuaternion(directionQuaternion, directionCosines);
////	directionTest
////	lpQuat.slerp(directionQuat4,alpha)
////	directionQuat4.
//}
//Vector3d svdConstellation(Matrix3Xd& constellation)
//{
//	Matrix3d constellationMatrix;
//	constellationMatrix.setZero();
//	for (int i = 0; i < constellation.cols(); i++)
//		constellationMatrix += constellation.col(i) * constellation.col(i).transpose();
//	JacobiSVD<Matrix3d> svdy(constellationMatrix);
//	return Vector3d(svdy.singularValues());
//}
void MovingPoseEstimator::setConstellation(Matrix3Xd constellation, double priorPositionWeight,
											double priorOrientationWeight, double timestep,
											double velocityPriorWeight)
{
	numPoints = constellation.cols();
	regressor.setZero(3 * numPoints + 15, 15);
	observations.setZero(3 * numPoints + 15);
	linearSolutionMatrix.setZero(15, 3 * numPoints + 15);
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
	regressor.block(numPoints * 3, 0, 15, 15) << Matrix15d::Identity();
	regressor.block(numPoints * 3, 12, 3, 3) << Matrix3d::Identity() * timestep;
	VectorXd weights(numPoints * 3 + 15);
	weights.setOnes();
	weights.segment(3 * numPoints, 3) *= priorPositionWeight;
	weights.segment(3 * numPoints + 3, 9) *= priorOrientationWeight;
	weights.segment(3 * numPoints + 12, 3) *= velocityPriorWeight;
	this->symmetricRegressionProduct = regressor.transpose() * weights.asDiagonal() * regressor;
	this->symmetricRegressionProductInverse = symmetricRegressionProduct.inverse();
	linearSolutionMatrix << symmetricRegressionProductInverse * regressor.transpose() * weights.asDiagonal();
//	std::cout << "regressor with priors\n" << regressor << std::endl;
//	std::cout << "solutionMatrix with priors\n" << linearSolutionMatrix << std::endl;
}

void MovingPoseEstimator::solve(Matrix3Xd observedPoints, Vector3d& expectedPosition, Matrix3d& expectedOrientation,
								Vector3d& expectedVelocity)
{
	unsigned numPoints = observedPoints.cols();
	populateVectorXdFromMatrix3Xd(observations, observedPoints);
	observations.segment(3 * numPoints, 3) << expectedPosition;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			observations(3 * numPoints + 3 + i * 3 + j) = expectedOrientation(i, j);
	observations.tail(3) << Vector3d::Zero();
	Vector15d result = linearSolutionMatrix * observations;
	center << result.segment(0, 3);
	rawDirectionCosines << result.segment(3, 3).transpose(), result.segment(6, 3).transpose(), result.segment(9, 3).transpose();
	velocityDeviation << result.tail(3);
	gu::vectorGeometry::closestQuaternionToMatrix(directionQuat4, rawDirectionCosines);
	gu::vectorGeometry::directionCosinesFromQuaternion(directionQuat4, directionCosines);
}

Quaterniond MovingPoseEstimator::getQuaternion()
{
	Quaterniond ret(directionQuat4);
	return ret;
}
Vector3d MovingPoseEstimator::getPoint()
{
	Vector3d ret(center);
	return ret;
}
Matrix3d MovingPoseEstimator::getDCM()
{
	Matrix3d ret(directionCosines);
	return ret;
}
//void populateVectorXdFromMatrix3Xd(VectorXd& vec, const Matrix3Xd& mat)
//{
//	for (int i = 0; i < mat.cols(); i++)
//		vec.segment(3 * i, 3) << mat.col(i);
//}
//void directionCosinesFromQuaternion(Vector4d quat, Matrix3d& dcm)
//{
//	quat.normalize();
//	double x = quat(0), y = quat(1), z = quat(2), w = quat(3);
//	dcm << 1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w), //
//	2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w), //
//	2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y);
//
//}
//Vector4d versor(Vector3d rot, double angle)
//{
//	rot.normalize();
//	double w = std::cos(angle / 2);
//	double p = std::sin(angle / 2);
//	Vector4d versor;
//	versor << p * rot(0), p * rot(1), p * rot(2), w;
//	return versor;
//}
//
//Vector3d unitVector(double x, double y, double z)
//{
//	Vector3d rot;
//	rot << x, y, z;
//	rot.normalize();
//	return rot;
//}

}
}
