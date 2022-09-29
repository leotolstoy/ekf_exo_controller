#ifndef TEST_NO_EIGEN__HPP
#define TEST_NO_EIGEN__HPP
#include "noEigenAttitudeEstimatorInterface.hpp"
#include "Daemon.hpp"
#include <iostream>
#include <unistd.h>
namespace gu
{
namespace attitudeEstimator
{
namespace noEigenInterface
{
void testNoEigen()
{
	HumePoseEstimator* poseEstimator = getLowPassEstimator();
	for (int i=0;i<1000;i++)
	{
		double x=poseEstimator->getPositionX(),y=poseEstimator->getPositionY(), z=poseEstimator->getPositionZ();
		double roll= poseEstimator->getRoll(), pitch=poseEstimator->getPitch(), yaw=poseEstimator->getYaw();
		std::cout<<"Position ["<<x<<", "<<y<<", "<<z<<"] and RPY {" << roll<<", "<<pitch<<", "<<yaw<<"}"<<std::endl;
		usleep(50000);
	}
}
}
}
}
#endif
