/*
 * AttitudeEstimator.hpp
 *
 *    Created on: May 13, 2014
 *        Author: Gray Thomas
 *   Affiliation: Human Centered Robotics Lab, University of Texas at Austin
 *
 */
#ifndef ATTITUDEESTIMATOR_HPP_
#define ATTITUDEESTIMATOR_HPP_
#include "Eigen/Dense"
#include "PhaseSpaceUDPReceiver.hpp"
#include <map>

using namespace Eigen;
namespace gu
{
namespace attitudeEstimator
{
typedef Matrix<double, Dynamic, 12> MatrixX12d;
typedef Matrix<double, 12, Dynamic> Matrix12Xd;
typedef Matrix<double, 12, 12> Matrix12d;
typedef Matrix<double, 3, Dynamic> Matrix3Xd;
typedef Matrix<double, 12, 1> Vector12d;
typedef Eigen::Quaterniond Quat4d;
class AttitudeEstimator
{
protected:
	Vector3d center;
	Matrix3d directionCosines;
	Matrix3d rawDirectionCosines;
//	Vector4d directionQuaternion;
	Quat4d directionQuat4;
	Matrix4d directionTest;
	MatrixX12d regressor;
	Matrix12d symmetricRegressionProduct;
	Matrix12d symmetricRegressionProductInverse;
	Matrix12Xd linearSolutionMatrix;
	VectorXd observations;
	Matrix3Xd constellation;
	int numPoints;

public:
	AttitudeEstimator();
	void solveDirectly(Matrix3Xd observedPoints);
	void solveDirectlyWithPrior(const Matrix3Xd& observedPoints, const Vector3d& priorPosition,
								const Matrix3d& priorOrientation);
	void setConstellation(Matrix3Xd constellation);
	void setConstellationWithPrior(const Matrix3Xd& constellation, double priorPositionWeight,
									double priorOrientationWeight);
	Quaterniond getQuaternion();
	Vector3d getPoint();
	Matrix3d getDCM();
};

void sayHi(void);
void populateVectorXdFromMatrix3Xd(VectorXd& vec, const Matrix3Xd& mat);
void moveConstellationCentroidToOrigin(gu::attitudeEstimator::Matrix3Xd& constellation);
Vector3d svdConstellation(Matrix3Xd& constellation);
//Vector4d versor(Vector3d rot, double angle);
Matrix3Xd getDefaultConstellation();

}
}

#endif /* ATTITUDEESTIMATOR_HPP_ */
