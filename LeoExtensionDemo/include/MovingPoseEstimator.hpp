/*
 * MovingPoseEstimator.hpp
 *
 *    Created on: June 21, 2014
 *        Author: Gray Thomas
 *   Affiliation: Human Centered Robotics Lab, University of Texas at Austin
 *
 */
#ifndef MOVINGPOSEESTIMATOR_HPP_
#define MOVINGPOSEESTIMATOR_HPP_
#include "Eigen/Dense"
#include "PhaseSpaceUDPReceiver.hpp"
#include <map>

using namespace Eigen;
namespace gu
{
namespace attitudeEstimator
{
typedef Matrix<double, Dynamic, 15> MatrixX15d;
typedef Matrix<double, 15, Dynamic> Matrix15Xd;
typedef Matrix<double, 15, 15> Matrix15d;
typedef Matrix<double, 3, Dynamic> Matrix3Xd;
typedef Matrix<double, 15, 1> Vector15d;
typedef Eigen::Quaterniond Quat4d;
class MovingPoseEstimator
{
protected:
	Vector3d center;
	Vector3d velocityDeviation;
	Matrix3d directionCosines;
	Matrix3d rawDirectionCosines;
	Vector4d directionQuaternion;
	Quat4d directionQuat4;
	Matrix4d directionTest;
	MatrixX15d regressor;
	Matrix15d symmetricRegressionProduct;
	Matrix15d symmetricRegressionProductInverse;
	Matrix15Xd linearSolutionMatrix;
	VectorXd observations;
	Matrix3Xd constellation;
	int numPoints;

public:
	MovingPoseEstimator();
	void solve(Matrix3Xd observedPoints, Vector3d& expectedPosition, Matrix3d& expectedOrientation,
				Vector3d& expectedVelocity);
	void setConstellation(Matrix3Xd constellation, double priorPositionWeight, double priorOrientationWeight,
							double timestep, double velocityPriorWeight);
	Quaterniond getQuaternion();
	Vector3d getPoint();
	Matrix3d getDCM();
};

//void populateVectorXdFromMatrix3Xd(VectorXd& vec, const Matrix3Xd& mat);
//void directionCosinesFromQuaternion(Vector4d quat, Matrix3d& dcm);
//Vector3d svdConstellation(Matrix3Xd& constellation);
//Vector4d versor(Vector3d rot, double angle);
//Vector3d unitVector(double x, double y, double z);
}
}

#endif /* ATTITUDEESTIMATOR_HPP_ */
