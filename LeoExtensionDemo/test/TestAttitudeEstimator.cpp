/*
 * TestAttitudeEstimator.cpp
 *
 *    Created on: May 13, 2014
 *        Author: Gray Thomas
 *   Affiliation: Human Centered Robotics Lab, University of Texas at Austin
 *
 */

#include "TestAttitudeEstimator.hpp"
#include "AttitudeEstimator.hpp"
#include "LowPassPositionFilter.hpp"
#include "HumePoseKalmanFilter.hpp"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "PhaseSpaceUDPReceiver.hpp"
#include "LowPassPositionFilter.hpp"
#include "HumePoseKalmanFilter.hpp"
#include "MovingPoseEstimator.hpp"
#include "LiveTest.hpp"
#include "TestWithVelocity.hpp"
#include <SpatialTransform.hpp>
#include <map>
#include <random>
#include <chrono>
#include "TestNoEigen.hpp"
using namespace Eigen;
using gu::attitudeEstimator::getNewTestConstellation1;
template<typename MatrixType>
void setMatrixGaussianRandom(MatrixType& input, double std, std::default_random_engine& gen)
{
	std::normal_distribution<double> normal(0, std);
	for (int r = 0; r < input.rows(); r++)
		for (int c = 0; c < input.cols(); c++)
			input(r, c) = normal(gen);
}
namespace gu{
namespace attitudeEstimator{

gu::attitudeEstimator::Matrix3Xd getTestConstellation1()
{
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << -574.67261, -577.24023, -598.6684, -735.65289, -717.95428, -433.29669, -449.78925, //
	1059.2797, 1350.1863, 1032.1157, 1211.9319, 1246.174, 1238.4542, 1215.3811, //
	2215.9297, 2052.9788, 1981.7324, 1977.489, 2159.1494, 2131.9507, 1962.0413;
	return ret;
}

gu::attitudeEstimator::Matrix3Xd getTestConstellation2()
{
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << -574.54407, -577.30151, -598.82025, -735.67188, -717.73248, -433.21027, -449.7926, //
	1059.3408, 1350.1976, 1032.052, 1211.918, 1246.1967, 1238.4634, 1215.3853, //
	2215.9617, 2053.0405, 1981.7603, 1977.4788, 2159.2549, 2131.8953, 1962.0439;
	return ret;
}
gu::attitudeEstimator::Matrix3Xd getTestConstellation3()
{
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << -574.6355, -577.22461, -598.63751, -735.67371, -717.96143, -433.28571, -449.83871, //
	1059.33, 1350.1753, 1032.1364, 1211.9207, 1246.1819, 1238.4315, 1215.3607, //
	2215.9614, 2052.9727, 1981.7128, 1977.4929, 2159.165, 2131.9299, 1962.0413;
	return ret;
}
gu::attitudeEstimator::Matrix3Xd getTestConstellation4()
{
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << -574.55157, -577.30817, -598.68915, -735.65839, -717.75049, -433.27271, -449.77878, //
	1059.3285, 1350.2037, 1032.1183, 1211.9233, 1246.225, 1238.4473, 1215.3892, //
	2215.9375, 2053.0317, 1981.7301, 1977.4828, 2159.2537, 2131.9304, 1962.0465;
	return ret;
}
gu::attitudeEstimator::Matrix3Xd getTestConstellation5()
{
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << -574.63708, -577.22656, -598.79803, -735.6803, -717.95422, -433.21259, -449.8313, //
	1059.3374, 1350.1722, 1032.0533, 1211.9052, 1246.1633, 1238.4587, 1215.3655, //
	2215.959, 2052.9771, 1981.7461, 1977.4965, 2159.1812, 2131.8884, 1962.057;
	return ret;
}
gu::attitudeEstimator::Matrix3Xd getNewTestConstellation1()
{
	MatrixXd retT(7, 3);
	retT << -585.411011, 1349.616577, 2144.938965, //
	-588.638977, 1615.051514, 1959.498047, //
	-605.916382, 1287.638550, 1893.115479, //
	-739.386719, 1470.277588, 1894.559692, //
	-716.418762, 1505.585449, 2061.909424, //
	-442.641052, 1497.346924, 2037.680054, //
	-457.334229, 1474.048584, 1873.852905;
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << retT.transpose();
	return ret;
}

gu::attitudeEstimator::Matrix3Xd getNewTestConstellation2()
{
	gu::attitudeEstimator::Matrix3Xd ret(3, 7);
	ret << -585.406, -588.688, -605.96, -739.41, -716.418, -442.685, -457.343, 1349.59, 1615.1, 1287.7, 1470.33, 1505.58, 1497.4, 1474.11, 2144.9, 1959.51, 1893.12, 1894.59, 2061.88, 2037.7, 1873.86;
	return ret;
}
void directionCosinesFromQuaternion(Vector4d quat, Matrix3d& dcm)
{
	quat.normalize();
	double x = quat(0), y = quat(1), z = quat(2), w = quat(3);
	dcm << 1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w), //
	2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w), //
	2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y);

}
void hello(void)
{
	std::cout << "Hello from the attitude Estimator Library" << std::endl;
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	gu::attitudeEstimator::Matrix3Xd constellation;
	constellation.setZero(3, 8);
//	constellation << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
	int col = 0;
	for (double x = -1; x < 2; x += 2)
		for (double y = -1; y < 2; y += 2)
			for (double z = -1; z < 2; z += 2)
				constellation.col(col++) << x, y, z;

	attitudeEstimator.setConstellation(constellation);
	Matrix3d rot;
	Vector4d quat;
	quat << 0, 0, 0, 1; // identity quaternion
	quat << 1, 0, 0, 0; // 180 x rotation quaternion (negates y and z)
	quat << 0, 1, 0, 0; // 180 y rotation quaternion (negates x and z)
	quat << 0, 0, 1, 0; // 180 z rotation quaternion (negates x and y)
	quat << 0, 0, 2, 0; // should be normalized by the function call
	quat << 0.5, 0, 0, 0.5; // 90 x rotation quaternion (y=-z, z=y)
	quat << -0.5, 0, 0, 0.5; // -90 x rotation quaternion (y=z, z=-y)
	quat << 0.1, 0, 0, 1.0; // slight x rotation quaternion
	Quaterniond versor1 = gu::vectorGeometry::versor(gu::vectorGeometry::unitVector(0.2, 0.4, -1.2),
			M_PI * 0.2);
	quat << versor1.x(), versor1.y(), versor1.z(), versor1.w();
	Vector3d vec;
	vec << 10.0, 0, 0;
//	vec.asDiagonal()* gu::attitudeEstimator::Matrix3Xd::Ones()

	std::cout << "\n" << quat << std::endl;
	directionCosinesFromQuaternion(quat, rot);

	gu::attitudeEstimator::Matrix3Xd noiseMatrix(3, 8);
	std::default_random_engine gen(234224L);
	setMatrixGaussianRandom(noiseMatrix, 0.1, gen);
	gu::attitudeEstimator::Matrix3Xd points = rot * constellation
			+ vec.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 8) + noiseMatrix;
	std::cout << "\n" << constellation << "\n" << std::endl;
	std::cout << points << std::endl;

	attitudeEstimator.solveDirectly(points);
	std::cout << "\nquat=\n" << quat << std::endl;

}
void testWithRealData(void)
{
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	gu::attitudeEstimator::Matrix3Xd constellation = getTestConstellation1();
	Vector3d originShift = constellation.rowwise().sum() * (1.0 / 7.0);
	constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);

	attitudeEstimator.setConstellation(constellation);
	attitudeEstimator.solveDirectly(
			getTestConstellation2() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7));
	attitudeEstimator.solveDirectly(
			getTestConstellation3() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7));
	attitudeEstimator.solveDirectly(
			getTestConstellation4() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7));
}

void testWithRealDataAndPriors(void)
{
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	gu::attitudeEstimator::Matrix3Xd constellation = getTestConstellation1();
	Vector3d originShift = constellation.rowwise().sum() * (1.0 / 7.0);
	constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);

	attitudeEstimator.setConstellation(constellation);
	attitudeEstimator.solveDirectly(
			getTestConstellation2() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7));
	Matrix3d DCM(attitudeEstimator.getDCM());
	Vector3d position(attitudeEstimator.getPoint());
	attitudeEstimator.setConstellationWithPrior(constellation, 7, 82000.0);
	attitudeEstimator.solveDirectlyWithPrior(
			getTestConstellation3() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7),
			position, DCM);
//	attitudeEstimator.solveDirectly(
//			getTestConstellation4() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7));
}

void testWithNewRealData(void)
{
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	gu::attitudeEstimator::Matrix3Xd constellation = getNewTestConstellation1();
	Vector3d originShift = constellation.rowwise().sum() * (1.0 / 7.0);
	constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);

	attitudeEstimator.setConstellation(constellation);
	attitudeEstimator.solveDirectly(
			getNewTestConstellation2() - originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7));
}
void testYawPitchRollConverter(void)
{
	std::cout << "testYawPitchRollConverter:" << std::endl;
	double yaw = -M_PI; //PUNT: setup actual asserts within some real, cmake compatible testing framework. ie.
	//# For make-based builds, defines make target named test.
	//# For Visual Studio builds, defines Visual Studio project named RUN_TESTS.
	//enable_testing()
	//Compile an executable that will run your unit tests and link it with gtest and gtest_main:
	//
	//add_executable(runUnitTests
	//    project1_unittests.cpp
	//)
	//target_link_libraries(runUnitTests gtest gtest_main)
	//Add a test which runs this executable:
	//
	//add_test(
	//    NAME runUnitTests
	//    COMMAND runUnitTests
	//)
	double pitch = -1.3;
	double roll = 1.5;
	Eigen::Matrix3d mat1 = gu::vectorGeometry::matrixFromYPR(yaw, pitch, roll);
	Eigen::Matrix3d mat1Confirm =  gu::vectorGeometry::matrixFromYPRConfirm(yaw, pitch, roll);
	std::cout << "test mat1 = \n"<<mat1<<std::endl;
	std::cout << "conf mat1 = \n"<<mat1Confirm<<std::endl;
	gu::vectorGeometry::packYawPitchRoll(yaw, pitch, roll, mat1);
	std::cout << "yaw = "<< yaw<<" pitch = "<<pitch<<" roll = "<<roll<<std::endl;

	std::cout << ":testYawPitchRollConverter" << std::endl;
}
}
}
int main(void)
{
//	hello();
//	testWithRealData();
//	testWithRealDataAndPriors();
//	testWithNewRealData();
//	testYawPitchRollConverter();
	testWithVelocity();
//	testLive();
//	gu::attitudeEstimator::noEigenInterface::testNoEigen();
//	pause();

	std::cout << "I guess it worked then" << std::endl;
	return 0;
}
