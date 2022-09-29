#ifndef DELAY_KALMAN_SMOOTHING_IMU_PHASESPACE_FILTER__HPP
#define DELAY_KALMAN_SMOOTHING_IMU_PHASESPACE_FILTER__HPP
#include <Eigen/Dense>
//#include <PhaseSpaceUDPReceiver.hpp>
#include <AttitudeEstimator.hpp>
#include <Frame.hpp>
#include <FrameVector.hpp>
#include "noEigenAttitudeEstimatorInterface.hpp"
#include <SpatialTransform.hpp>
#include "MatrixListener.hpp"
#include "MatrixReceiver.hpp"
#include <map>
#include <list>

namespace gu
{
namespace attitudeEstimator
{
// Pseudo code for update:
// get new LED data:
// 1: find the delayth previous poseHistoryElement in poseHistory
// 2: calculate the expected transform from camera to constellation, using the transform from this history element
// 2 and haif: calculate the expected constellation velocity, using the velocity from this history element
// 3: solve for the best fit camera to constellation pose, using heavy regularization to this good guess (steady state Kalman filter)
// 3.5: divide the update in position between a position update and a velocity update.
// 4: check result for crazy stuff, and ignore the update if it would introduce an unreasonable value (which should never happen)
// 5: this result is used to compute a new bodyPoseLastLED and bodyVelocityLastLED
// 6: remove imuData before this time, and clear poseHistory
// 7: recalculate the poseHistory using bodyPoseLastLED, bodyVelocityLastLED by iterating over the imuData.

// ongoing: as more imu data comes in, keep the poseHistory up to date
// when requested: return the end transform in poseHistory when asked for the pose.
// when velocity is requested, return the end velocity in poseHistory
// when angular velocity is requested, return the angular velocity from the most recent imu, converted to the appropriate frame
// when acceleration is requested, return the most recent imu entry
// when angular acceleration is requested, we will use a simple filter to obtain the derivative data from the imu.
extern char bodyFrameName[];
extern char imuFrameName[];
extern char cameraFrameName[];
extern char ledFrameName[];

using gu::vectorGeometry::worldFrame;
using gu::vectorGeometry::Frame;
using gu::vectorGeometry::FramePoint;
using gu::vectorGeometry::FrameVector;
using gu::vectorGeometry::SpatialTransform;

typedef Frame<bodyFrameName> bodyFrame;
typedef Frame<imuFrameName> imuFrame;
typedef Frame<ledFrameName> ledFrame;
typedef Frame<cameraFrameName> cameraFrame;

char const* greet();

struct IMUDataPoint
{
	FrameVector<imuFrame> acceleration;
	FrameVector<imuFrame> angularVelocity;
	IMUDataPoint(const FrameVector<imuFrame>& acceleration, const FrameVector<imuFrame>& angularVelocity) :
				acceleration(acceleration),
				angularVelocity(angularVelocity)
	{
	}
};

struct PoseHistoryElement
{
	SpatialTransform<worldFrame, bodyFrame> pose;
	FrameVector<worldFrame> velocity;
	PoseHistoryElement() :
				pose(),
				velocity()
	{
	}
};
class HumePoseKalmanFilter: public gu::humeEstimation::MatrixListener<7>, public noEigenInterface::HumePoseEstimator
{
	int count = 0;
	std::map<uint32_t, gu::attitudeEstimator::AttitudeEstimator> alternateEstimators;
	std::map<uint32_t, Vector3d> constellationSvds;
	Eigen::Matrix<double, 7, 7> selector;
	Eigen::Vector3d currentPosition;
	Eigen::Matrix3d currentDirectionCosineMatrix;
	Eigen::Quaterniond currentQuaternion;
	double priorPointWeight = 7.0, priorMatrixWeight = 80000.0;
	unsigned numberOfPoints = 0;
	double averageVisibleLEDs = 7;
	Eigen::Vector3d averageSingularValuesOfConstellation;
	std::vector<gu::humeEstimation::MatrixReceiver<7>*> managedReceivers;

private:
	//data
	FrameVector<worldFrame> inertialFrameAcceleration;
	int delay;
	double timestep;
	bool ledUpdating;
private:
	//variables
	double x, y, z, roll, pitch, yaw;
	SpatialTransform<worldFrame, bodyFrame> bodyPoseLastLED;
	FrameVector<worldFrame> bodyVelocityLastLED;
	std::list<IMUDataPoint> imuData;
	std::list<PoseHistoryElement> poseHistory;
	void setupPoses();
//	void printMessage();
	void updateSelector(const Matrix<bool, 1, 7> &valid_leds);
	void updateSVDMetaData(int32_t validBits);
	SpatialTransform<worldFrame, ledFrame> solveLEDSubproblem(
			const SpatialTransform<worldFrame, ledFrame>& estimatedHistoricLEDPoseInWorldFrame,
			const Matrix<double, 3, 7>& led_mesaurements, const Matrix<bool, 1, 7>& led_selector);
public:
	SpatialTransform<bodyFrame, ledFrame> ledPose;
	SpatialTransform<bodyFrame, imuFrame> imuPose;
	SpatialTransform<worldFrame, cameraFrame> cameraPose;
	Matrix3Xd ledsInLEDFrame;
	// convenient methods
	HumePoseKalmanFilter(double priorPointWeight, double priorMatrixWeight, const Matrix3d& initialPose,
							const Vector3d& initialPosition, const Vector3d& initialVelocity);
	~HumePoseKalmanFilter();
	void handleUpdate(const Matrix<double, 3, 7>&, const Matrix<bool, 1, 7>&);
	inline void add_dependent_receiver(gu::humeEstimation::MatrixReceiver<7>* receiver)
	{
		managedReceivers.push_back(receiver);
	}
	const SpatialTransform<worldFrame, bodyFrame>& getTransform() const
	{
		return this->poseHistory.back().pose;
	}
	const FrameVector<worldFrame>& getVelocity() const
	{
		return bodyVelocityLastLED;
	}
	void setLEDPose(const vectorGeometry::SpatialTransform<bodyFrame, ledFrame>& transform)
	{
		ledPose.setPosition(transform.getPosition());
		ledPose.setOrientation(transform.getQuaternion());
	}
	void setIMUPose(const vectorGeometry::SpatialTransform<bodyFrame, imuFrame>& transform)
	{
		imuPose = transform;
	}
//	void spoofLEDData(const Matrix3Xd& constellation, Matrix<bool, 1, 7> selector);

protected:
	void lazyInitializeAttitudeEstimator();
	void handleUpdate();
	void integrateIMU(const IMUDataPoint& point);
public:
	// no eigen interface methods
	void updateIMU(double accelX, double accelY, double accelZ, double omegaX, double omegaY, double omegaZ);

	double getPositionX()
	{
		return x;
	}
	double getPositionY()
	{
		return y;
	}
	double getPositionZ()
	{
		return z;
	}
	double getRoll()
	{
		return roll;
	}
	double getPitch()
	{
		return pitch;
	}
	double getYaw()
	{
		return yaw;
	}
	const Eigen::Quaterniond& getOrientation() const
	{
		return currentQuaternion;
	}
	const Eigen::Matrix3d& getDCM() const
	{
		return currentDirectionCosineMatrix;
	}
	double getAverageNumberOfLEDs() const
	{
		return averageVisibleLEDs;
	}
	Eigen::Vector3d getAverageSVDForOrientationEstimation() const
	{
		return Eigen::Vector3d(averageSingularValuesOfConstellation);
	}

};
}
}
#endif
