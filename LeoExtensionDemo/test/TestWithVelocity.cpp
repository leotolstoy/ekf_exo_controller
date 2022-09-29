#include "TestWithVelocity.hpp"
#include "TestAttitudeEstimator.hpp"
#include <SpatialTransform.hpp>
#include "HumePoseKalmanFilter.hpp"
#include "Frame.hpp"
#include <iostream>
#include <random>
using gu::attitudeEstimator::getNewTestConstellation1;
using gu::attitudeEstimator::bodyFrame;
using gu::vectorGeometry::worldFrame;
using gu::attitudeEstimator::cameraFrame;
using gu::attitudeEstimator::imuFrame;
using gu::attitudeEstimator::ledFrame;
using gu::vectorGeometry::SpatialTransform;
using gu::vectorGeometry::FrameVector;
gu::attitudeEstimator::Matrix3Xd corrupt(gu::attitudeEstimator::Matrix3Xd& input, double noise, long seed)
{
	std::default_random_engine engine(seed + 123432);
	std::normal_distribution<double> noiseDist(0.0, noise);
	gu::attitudeEstimator::Matrix3Xd output(input);
	for (int r = 0; r < input.rows(); r++)
		for (int c = 0; c < input.cols(); c++)
			output(r, c) += noiseDist(engine);
	return output;
}

gu::attitudeEstimator::Matrix3Xd move(const gu::attitudeEstimator::Matrix3Xd& input, Vector3d move)
{
	gu::attitudeEstimator::Matrix3Xd output(input);
	output.colwise() += move;
	return output;
}
gu::attitudeEstimator::Matrix3Xd transform(const gu::attitudeEstimator::Matrix3Xd& input, const Vector3d& move, const Matrix3d& dcm);
gu::attitudeEstimator::Matrix3Xd transform(const gu::attitudeEstimator::Matrix3Xd& input, const Vector3d& move, const Matrix3d& dcm)
{
	gu::attitudeEstimator::Matrix3Xd output = dcm * input;
	output.colwise() += move;
	return output;
}

void testOneUpdateWithSimpleNoise(void)
{
	gu::attitudeEstimator::Matrix3Xd startingConstellation = getNewTestConstellation1();
	gu::attitudeEstimator::moveConstellationCentroidToOrigin(startingConstellation);
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	attitudeEstimator.setConstellationWithPrior(startingConstellation, 7, 80000);
	Vector3d positionGuess;
	Matrix3d dcmGuess;
	dcmGuess << Matrix3d::Identity();
	for (int i = 1; i < 10; i++)
	{
		attitudeEstimator.solveDirectlyWithPrior(corrupt(startingConstellation, 2.0, i), positionGuess, dcmGuess);
//		std::cout << attitudeEstimator.getQuaternion();
//		std::cout << std::endl;
		double yaw, pitch, roll;
		gu::vectorGeometry::packYawPitchRoll(yaw, pitch, roll, attitudeEstimator.getDCM());
		std::cout << "yaw = " << yaw * 180 / M_PI << "\u00b0" " pitch = " << pitch * 180 / M_PI << "\u00b0" " roll = " << roll * 180 / M_PI
				<< "\u00b0" << std::endl;
	}
}
void testOneUpdateWithOffset(void)
{
	std::cout << "now with an offset" << std::endl;
	gu::attitudeEstimator::Matrix3Xd startingConstellation = getNewTestConstellation1();
	gu::attitudeEstimator::moveConstellationCentroidToOrigin(startingConstellation);
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	attitudeEstimator.setConstellationWithPrior(startingConstellation, 7, 80000);
	Vector3d actualPosition;
	actualPosition << 100.0, 0, 0;
	Vector3d positionGuess;
	positionGuess << 90.0, 0, 0;
	Matrix3d dcmGuess;
	gu::attitudeEstimator::Matrix3Xd movedConstellation = move(startingConstellation, actualPosition);
	dcmGuess << Matrix3d::Identity();
	for (int i = 1; i < 10; i++)
	{
		attitudeEstimator.solveDirectlyWithPrior(corrupt(movedConstellation, 2.0, i), positionGuess, dcmGuess);
//		std::cout << attitudeEstimator.getQuaternion();
//		std::cout << std::endl;
		double yaw, pitch, roll;
		gu::vectorGeometry::packYawPitchRoll(yaw, pitch, roll, attitudeEstimator.getDCM());
		Vector3d point = attitudeEstimator.getPoint();
		std::cout << "yaw = " << yaw * 180 / M_PI << "\u00b0" " pitch = " << pitch * 180 / M_PI << "\u00b0" " roll = " << roll * 180 / M_PI
				<< "\u00b0 point = " << point.transpose();
		std::cout << std::endl;
	}
}
inline double deg(double deg)
{
	return deg * M_PI / 180.0;
}
void testOneUpdateWithRotation(void)
{
	std::cout << "now with a transform" << std::endl;
	gu::attitudeEstimator::Matrix3Xd startingConstellation = getNewTestConstellation1();
	gu::attitudeEstimator::moveConstellationCentroidToOrigin(startingConstellation);
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;
	attitudeEstimator.setConstellationWithPrior(startingConstellation, 7, 80000);
	Vector3d actualPosition;
	actualPosition << 100.0, 0, 0;
	Vector3d positionGuess;
	positionGuess << 90.0, 0, 0;
	Matrix3d actualDCM = gu::vectorGeometry::matrixFromYPR(deg(10.0), deg(10.0), deg(10.0));
	Matrix3d dcmGuess = gu::vectorGeometry::matrixFromYPR(deg(10.0), deg(10.0), deg(10.0));
	gu::attitudeEstimator::Matrix3Xd movedConstellation = transform(startingConstellation, actualPosition, actualDCM);
//	dcmGuess << Matrix3d::Identity();
	for (int i = 1; i < 10; i++)
	{
		attitudeEstimator.solveDirectlyWithPrior(corrupt(movedConstellation, 2.0, i), positionGuess, dcmGuess);
//		std::cout << attitudeEstimator.getQuaternion();
//		std::cout << std::endl;
		double yaw, pitch, roll;
		gu::vectorGeometry::packYawPitchRoll(yaw, pitch, roll, attitudeEstimator.getDCM());
		Vector3d point = attitudeEstimator.getPoint();
		std::cout << "yaw = " << yaw * 180 / M_PI << "\u00b0" " pitch = " << pitch * 180 / M_PI << "\u00b0" " roll = " << roll * 180 / M_PI
				<< "\u00b0 point = " << point.transpose();
		std::cout << std::endl;
	}
}

void testOneUpdateWithExtraTransforms(void)
{
	std::cout << "now with extra transforms" << std::endl;

	Vector3d cameraPosition, baseOriginInCameraFrame;
	Matrix3d cameraDCM;
	cameraDCM << 0, 0, -1, -1, 0, 0, 0, 1, 0;
	baseOriginInCameraFrame << -528, 0, 3364;
	cameraPosition = -cameraDCM * baseOriginInCameraFrame;

	Vector3d imuPositionInLedFrame;
	imuPositionInLedFrame << 20, 50, 100;
	Matrix3d imuDcmInLedFrame;
	imuDcmInLedFrame << 0, -1, 0, -1, 0, 0, 0, 0, -1;

	gu::attitudeEstimator::Matrix3Xd startingConstellation = transform(getNewTestConstellation1(), cameraPosition, cameraDCM);
	gu::attitudeEstimator::moveConstellationCentroidToOrigin(startingConstellation);
	gu::attitudeEstimator::AttitudeEstimator attitudeEstimator;

	attitudeEstimator.setConstellationWithPrior(startingConstellation, 7, 80000);

	Vector3d actualPosition;
	actualPosition << 100.0, 0, 0;
	Vector3d positionGuess;
	positionGuess << 90.0, 0, 0;
	Matrix3d actualDCM = gu::vectorGeometry::matrixFromYPR(deg(10.0), deg(10.0), deg(10.0));
	Matrix3d dcmGuess = gu::vectorGeometry::matrixFromYPR(deg(10.0), deg(10.0), deg(10.0));
	gu::attitudeEstimator::Matrix3Xd movedConstellation = transform(startingConstellation, actualPosition, actualDCM);
//	dcmGuess << Matrix3d::Identity();
	for (int i = 1; i < 10; i++)
	{
		attitudeEstimator.solveDirectlyWithPrior(corrupt(movedConstellation, 2.0, i), positionGuess, dcmGuess);
//		std::cout << attitudeEstimator.getQuaternion();
//		std::cout << std::endl;
		double yaw, pitch, roll;
		gu::vectorGeometry::packYawPitchRoll(yaw, pitch, roll, attitudeEstimator.getDCM());
		Vector3d point = attitudeEstimator.getPoint();
		std::cout << "yaw = " << yaw * 180 / M_PI << "\u00b0" " pitch = " << pitch * 180 / M_PI << "\u00b0" " roll = " << roll * 180 / M_PI
				<< "\u00b0 point = " << point.transpose();
		std::cout << std::endl;
	}
}
template<typename BaseFrame, typename ChildFrame>
class FlawedRigidBodyTimeTrajectory
{
public:
	virtual ~FlawedRigidBodyTimeTrajectory(){}
	virtual SpatialTransform<BaseFrame,ChildFrame> getTransform(double t)=0;
	virtual FrameVector<BaseFrame> getV(double t)=0;
	virtual FrameVector<BaseFrame> getA(double t)=0;
	virtual FrameVector<BaseFrame> getOmega(double t)=0;
	virtual FrameVector<BaseFrame> getAlpha(double t)=0;
};
template<typename BaseFrame, typename ChildFrame>
class SinusoidalTrajectory: public FlawedRigidBodyTimeTrajectory<BaseFrame, ChildFrame>
{
	Vector3d direction;
	Vector3d axis;
	Vector3d initialPosition;
	double linearAmount;
	double rotationalAmount;
	double phaseShift;
	double frequency;
	double getPhase(double time)
	{
		return phaseShift + time * 2 * M_PI * frequency;
	}

	Quaterniond getRawTheta(double t)
	{
		return gu::vectorGeometry::versor(axis, rotationalAmount * sin(getPhase(t)));
	}
	Vector3d getRawOmega(double t)
	{
		Vector3d ret(axis);
		ret *= rotationalAmount * 2 * M_PI * frequency * cos(getPhase(t));
		return ret;
	}
	Vector3d getRawAlpha(double t)
	{
		Vector3d ret(axis);
		ret *= -rotationalAmount * 4 * M_PI * M_PI * frequency * frequency * sin(getPhase(t));
		return ret;
	}
	Vector3d getRawX(double t)
	{
		Vector3d ret(direction);
		ret *= linearAmount * sin(getPhase(t));
		return ret + initialPosition;
	}
	Vector3d getRawV(double t)
	{
		Vector3d ret(direction);
		ret *= linearAmount * cos(getPhase(t)) * 2 * M_PI * frequency;
		return ret;
	}
	Vector3d getRawA(double t)
	{
		Vector3d ret(direction);
		ret *= -linearAmount * sin(getPhase(t)) * 4 * M_PI * M_PI * frequency * frequency;
		return ret;
	}
public:
	explicit SinusoidalTrajectory( double x, double y, double z,
			double rx, double ry, double rz, double freq, double shift) :
					initialPosition()
	{
		direction << x, y, z;
		axis << rx, ry, rz;
		linearAmount = direction.norm();
		if (direction.norm() < 1e-14)
			direction << 1, 0, 0;
		direction.normalize();
		rotationalAmount = axis.norm();
		if (axis.norm() < 1e-14)
			axis << 1, 0, 0;
		axis.normalize();
		phaseShift = shift;
		frequency = freq;
		initialPosition.setZero();
	}
	virtual ~SinusoidalTrajectory(){}
	void setInitialPosition(const Eigen::Vector3d& vec)
	{
		initialPosition = vec;
	}
	virtual SpatialTransform<BaseFrame,ChildFrame> getTransform(double t)
	{
		SpatialTransform<BaseFrame,ChildFrame> transform;
		transform.setOrientation(getRawTheta(t));
		transform.setPosition(getRawX(t));
		return transform;
	}
	virtual FrameVector<BaseFrame> getV(double t)
	{
		FrameVector<BaseFrame> vectorV(getRawV(t));
		return vectorV;
	}
	virtual FrameVector<BaseFrame> getA(double t)
	{
		FrameVector<BaseFrame> vectorA(getRawA(t));
		return vectorA;
	}
	virtual FrameVector<BaseFrame> getOmega(double t)
	{
		FrameVector<BaseFrame> vectorOmega(getRawOmega(t));
		return vectorOmega;
	}
	virtual FrameVector<BaseFrame> getAlpha(double t)
	{
		FrameVector<BaseFrame> vectorAlpha(getRawAlpha(t));
		return vectorAlpha;
	}

};

//PUNT: quaternion slerp bezier curves
//PUNT: quaternion differentiation in vectorGeometry library
using namespace gu::vectorGeometry;
void testHumePoseKalmanFilter(void)
{
	std::cout << "testHumePoseKalmanFilter:" << std::endl;

	SpatialTransform<worldFrame,cameraFrame> cameraPose;
	SpatialTransform<worldFrame,bodyFrame> bodyPose;
	SpatialTransform<bodyFrame,ledFrame> ledPose;
	SpatialTransform<bodyFrame,imuFrame> imuPose;
	// Parentage tree for frames is currently meaningless. Parentage tree for transforms is critical.

	Vector3d cameraPosition, baseOriginInCameraFrame;
	Matrix3d cameraDCM;
	cameraDCM << 0, 0, -1, -1, 0, 0, 0, 1, 0;
	baseOriginInCameraFrame << -.528, 0, 3.364;
	// This is just a guess at the moment.
	cameraPosition = -cameraDCM * baseOriginInCameraFrame;
	cameraPose.setOrientation(cameraDCM);
	cameraPose.setPosition(cameraPosition);

	gu::attitudeEstimator::Matrix3Xd ledsInCameraFrame = 1e-3 * getNewTestConstellation1();
	gu::attitudeEstimator::Matrix3Xd ledsInWorldFrame = cameraPose * ledsInCameraFrame;
	Vector3d constellationCentroidInCamera = ledsInCameraFrame.rowwise().sum() * (1.0 / 7.0);
	Vector3d constellationCentroidInWorld = ledsInWorldFrame.rowwise().sum() * (1.0 / 7.0);
	SpatialTransform<worldFrame,ledFrame> initialLEDPose;
	initialLEDPose.setPosition(constellationCentroidInWorld); // initial led frame is at the centroid of the leds
	initialLEDPose.setOrientation(Matrix3d::Identity()); // initial led frame is aligned with world axes
	Matrix3Xd ledsInLEDFrame = initialLEDPose.inverse() * ledsInWorldFrame;

	std::cout << "\nAssertion: equal \n\t(";
//	std::cout << ((Vector3d) (cameraPose * constellationCentroidInCamera)).transpose();
	std::cout << ") = cameraPose*constellationCentroidInCamera \n\t(";
	std::cout << constellationCentroidInWorld.transpose();
	std::cout << ") = constellationCentroidInWorld ";
	std::cout << std::endl;

	std::cout << "\nNew constellation: leds in led frame:\n" << ledsInLEDFrame << std::endl;

	imuPose.setOrientation(0, -1, 0, -1, 0, 0, 0, 0, -1);
	imuPose.setPosition(.020, .050, 0.100);

	ledPose.setOrientation(1.0, 0.0, 0.0, 0.0);
	ledPose.setPosition(0.0, 0.0, 0.0); //Body frame is coincident with led frame for now

	FrameVector<worldFrame> gravitationalAcceleration( 0.0, 0.0, -9.81);

	SinusoidalTrajectory<worldFrame,imuFrame> imuTrajectory(2.0, 2.0, 0.0, 0.2, 0.3, 0.4, 0.1, 0);
	imuTrajectory.setInitialPosition(initialLEDPose.getPosition());
	//generate rotation error:
	Matrix3d startingError = matrixFromYPR(0.5, 0.0, 0.0);
	gu::attitudeEstimator::HumePoseKalmanFilter filter(0.01, 0.080,
			(imuTrajectory.getTransform(0) * imuPose.inverse()).getDCM() * startingError,
			(imuTrajectory.getTransform(0) * imuPose.inverse()).getPosition(), imuTrajectory.getV(0).getVector());

	double startTime = 0;
	double endTime = 0.05;
	double IMUFreq = 1000;
	double LEDFreq = 460;
	double ledCorruption = 0.01;

	std::cout << "Generating data" << std::endl;
	int numIMUData = (endTime - startTime) * IMUFreq;

	for (int i = 0; i < numIMUData; i++)
	{
		double t = startTime + 1.0 / IMUFreq * i;
		if (i % 2 == 0) //TODO: try non-synchronous update
		{ //TODO: try imu first initialization
			std::cout << "============================ LEDGEN at t=" << t << " ======\n";
			//LED update too.
			bodyPose = imuTrajectory.getTransform(t - 0.015) * imuPose.inverse(); //TODO: identify real delay
			Matrix3d ledFrameIntoCameraFrame;
			Vector3d ledFrameOriginInCameraFrame;
			SpatialTransform<cameraFrame,ledFrame> ledFrameInCameraFrame = cameraPose.inverse() * bodyPose * ledPose;
			gu::attitudeEstimator::Matrix3Xd observedLEDSinCamreraFrame = 1000.0
					* (cameraPose.inverse() * bodyPose * ledPose * ledsInLEDFrame);
			Matrix3Xd rawLEDData = corrupt(observedLEDSinCamreraFrame, ledCorruption, i + 100); //convert to mm;
			std::cout << "\trawLEDData[" << i << "]=\n" << rawLEDData << std::endl;
			Matrix<bool, 1, 7> bools;
			bools << true, true, true, true, true, true, true;
			Matrix<double,3,7> data(rawLEDData);
			filter.handleUpdate(data, bools);
			std::cout << filter.getTransform() << " = filter.getTransform()" << std::endl;
			std::cout << imuTrajectory.getTransform(t) * imuPose.inverse() << " = imuTrajectory.getTransform(t) * imuPose.inverse()"
					<< std::endl;

		}
		std::cout << "============================ IMUGEN at t=" << t << " ======\n";
		bodyPose = imuTrajectory.getTransform(t) * imuPose.inverse();
		FrameVector<imuFrame> imuA = imuPose.inverse() * bodyPose.inverse() * (imuTrajectory.getA(t) - gravitationalAcceleration);
		FrameVector<imuFrame> imuOmega = imuPose.inverse() * bodyPose.inverse() * (imuTrajectory.getOmega(t));
		Vector3d imuAccel = imuA.getVector() * 1.0 / 9.81;	// unit conversion
		Vector3d imuAngVel = imuOmega.getVector();
		std::cout << "\trawIMUData : omega=<" << imuAngVel.transpose() << ">, accel=<" << imuAccel.transpose() << ">\n";
		filter.updateIMU((double) imuAccel[0], (double) imuAccel[1], (double) imuAccel[2], (double) imuAngVel[0], (double) imuAngVel[1],
				(double) imuAngVel[2]);
		std::cout << filter.getTransform() << " = filter.getTransform()" << std::endl;
		std::cout << imuTrajectory.getTransform(t) * imuPose.inverse() << " = imuTrajectory.getTransform(t) * imuPose.inverse()"
				<< std::endl;
	}

}

void testWithVelocity(void)
{
	std::cout << "testWithVelocity:" << std::endl;
//	testOneUpdateWithSimpleNoise();
//	testOneUpdateWithOffset();
//	testOneUpdateWithRotation();
//	testOneUpdateWithExtraTransforms();
	testHumePoseKalmanFilter();

	std::cout << ":testWithVelocity" << std::endl;
}
