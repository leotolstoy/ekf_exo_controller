#include "HumePoseKalmanFilter.hpp"
#include "AttitudeEstimator.hpp"
#include "HumePoseKalmanFilter.hpp"
#include "SpatialTransform.hpp"

namespace gu
{
namespace attitudeEstimator
{
char worldFrameName[] = "World Frame";
char bodyFrameName[] = "Body Frame";
char imuFrameName[] = "IMU Frame";
char cameraFrameName[] = "Camera Frame";
using namespace gu::vectorGeometry;

noEigenInterface::HumePoseEstimator* noEigenInterface::getDefaultPoseEstimator()
{
	// better values: 25.0, 0.4;
	HumePoseKalmanFilter* filter = new HumePoseKalmanFilter(0.01, 0.080, Matrix3d::Identity(), Vector3d::Zero(), Vector3d::Zero());
	gu::humeEstimation::MatrixReceiver<7> *receiver= new gu::humeEstimation::MatrixReceiver<7>();
	receiver->add_listener(filter);
	filter->add_dependent_receiver(receiver);// Set proper destructor for the receiver.
	receiver->start();
	return filter;
}

char const* greet()
{
	return "hello world from HumePoseKalmanFilter.cpp.";
}

//uta::attitudeEstimator::Matrix3Xd getDefaultConstellation()
//{
//	MatrixXd retT(7, 3);
//	retT << -585.411011, 1349.616577, 2144.938965, //
//	-588.638977, 1615.051514, 1959.498047, //
//	-605.916382, 1287.638550, 1893.115479, //
//	-739.386719, 1470.277588, 1894.559692, //
//	-716.418762, 1505.585449, 2061.909424, //
//	-442.641052, 1497.346924, 2037.680054, //
//	-457.334229, 1474.048584, 1873.852905;
//	uta::attitudeEstimator::Matrix3Xd ret(3, 7);
//	ret << retT.transpose();
//	return ret;
//}
int getDefaultDelay()
{
	return 4; // unestimated
}

HumePoseKalmanFilter::HumePoseKalmanFilter(double priorPointWeight, double priorMatrixWeight, const Matrix3d& initialOrientation,
		const Vector3d& initialPosition, const Vector3d& initialVelocity) :
				alternateEstimators(),
				selector(),
				currentPosition(),
				currentDirectionCosineMatrix(),
				currentQuaternion(),
				priorPointWeight(priorPointWeight),
				priorMatrixWeight(priorMatrixWeight),
				numberOfPoints(0),
				averageVisibleLEDs(7),
				averageSingularValuesOfConstellation(),
				cameraPose(),
				ledPose(),
				imuPose(),
				ledsInLEDFrame(),
				inertialFrameAcceleration(0.0, 0.0, -9.8),
				delay(2),
				timestep(0.001),
				ledUpdating(false),
				bodyPoseLastLED(),
				bodyVelocityLastLED(),
				imuData(),
				poseHistory(),
				x(0),
				y(0),
				z(0),
				roll(0),
				pitch(0),
				yaw(0),
				managedReceivers()
{
	this->priorPointWeight = priorPointWeight;
	this->priorMatrixWeight = priorMatrixWeight;
	// setup poses

	Vector3d cameraPosition, baseOriginInCameraFrame;
	Matrix3d cameraDCM;
	cameraDCM << 0, 0, -1, -1, 0, 0, 0, 1, 0;
	baseOriginInCameraFrame << -.528, 0, 3.364;
	// This is just a guess at the moment.
	cameraPosition = -cameraDCM * baseOriginInCameraFrame;
	cameraPose.setOrientation(cameraDCM);
	cameraPose.setPosition(cameraPosition);

	gu::attitudeEstimator::Matrix3Xd ledsInCameraFrame = 1e-3 * getDefaultConstellation();
	gu::attitudeEstimator::Matrix3Xd ledsInWorldFrame = cameraPose * ledsInCameraFrame;
	Vector3d constellationCentroidInCamera = ledsInCameraFrame.rowwise().sum() * (1.0 / 7.0);
	Vector3d constellationCentroidInWorld = ledsInWorldFrame.rowwise().sum() * (1.0 / 7.0);
	SpatialTransform<worldFrame, ledFrame> initialLEDPose;
	initialLEDPose.setPosition(constellationCentroidInWorld); // initial led frame is at the centroid of the leds
	initialLEDPose.setOrientation(Matrix3d::Identity()); // initial led frame is aligned with world axes
	ledsInLEDFrame = initialLEDPose.inverse() * ledsInWorldFrame;
	imuPose.setOrientation(0, -1, 0, -1, 0, 0, 0, 0, -1);
	imuPose.setPosition(.020, .050, 0.100);
	ledPose.setOrientation(1.0, 0.0, 0.0, 0.0);
	ledPose.setPosition(0.0, 0.0, 0.0); //Body frame is coincident with led frame for now
	// \poses

	PoseHistoryElement startingPoseHistoryElement;
	startingPoseHistoryElement.pose.setPosition(initialPosition);
	startingPoseHistoryElement.pose.setOrientation(initialOrientation);
	startingPoseHistoryElement.velocity = FrameVector<worldFrame>(initialVelocity);
	poseHistory.push_back(startingPoseHistoryElement);

}
HumePoseKalmanFilter::~HumePoseKalmanFilter()
{
	for (auto it=managedReceivers.begin();it!=managedReceivers.end();++it)
	{
		delete (*it);
	}
}
void HumePoseKalmanFilter::handleUpdate(const Matrix<double,3,7>& leds,const Matrix<bool,1,7>& selector)
{
	ledUpdating = true;
	// find appropriate entry in poseHistory
	std::list<PoseHistoryElement>::iterator poseHistoryIterator = poseHistory.end();
	std::list<IMUDataPoint>::iterator imuDataIterator = imuData.end();
	for (int i = 0; i < delay; i++)
		if (--poseHistoryIterator == poseHistory.begin() || --imuDataIterator == imuData.begin())// are we sure these are always both executed?
			break;
	// the above returns the poseHistoryIterator and the imuDataIterator either to the beginning or back delay steps.
	// LED imu imu imu imu imu imu imu imu
	//  pos pos pos pos pos pos pos pos pos
	// poseHistory is one element larger than imuData, so we discard the imudata from the same time instant as the new LED data
	// as this imu point already effects the prior for that time instant
	imuData.erase(imuData.begin(), ++imuDataIterator); // I'm not sure this pointer/iterator stuff was really tested.

	SpatialTransform<worldFrame, ledFrame> estimatedHistoricLEDPoseInWorldFrame = poseHistoryIterator->pose * ledPose;
	SpatialTransform<worldFrame, ledFrame> newEstimatedHistoricLEDPoseInWorldFrame = solveLEDSubproblem(
			estimatedHistoricLEDPoseInWorldFrame, leds, selector);

	poseHistory.clear();
	//PUNT: pose history should mark the pose of the imu frame.
	PoseHistoryElement startingElement;
	startingElement.pose = newEstimatedHistoricLEDPoseInWorldFrame * ledPose.inverse();
	startingElement.velocity = FrameVector<worldFrame>(Vector3d::Zero());
	//PUNT: handle velocity properly
	poseHistory.push_back(startingElement);
	for (; imuDataIterator != imuData.end(); imuDataIterator++)
		integrateIMU(*imuDataIterator);
	//PUNT: pose history should mark the pose of the imu frame, so a transform will be necessary here.
	packYawPitchRoll(yaw, pitch, roll, poseHistory.back().pose.getDCM());
	Vector3d xyz(poseHistory.back().pose.getPosition());
	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
	ledUpdating = false;
}
void HumePoseKalmanFilter::handleUpdate()
{

}
//void HumePoseKalmanFilter::printMessage(void)
//{
//	message* msg = &mostRecentMessage;
//	Int32BitField* bitfield = (Int32BitField*) &mostRecentMessage.validBits;
//	bool validBools[7] =
//	{ (bool) bitfield->b0, (bool) bitfield->b1, (bool) bitfield->b2, (bool) bitfield->b3, (bool) bitfield->b4, (bool) bitfield->b5,
//			(bool) bitfield->b6 };
//	fprintf(stdout,
//			"Received Data: {%d,%d,%d, %d,%d,%d,%d} \n\t(%f %f %f)(%f %f %f)(%f %f %f)\n\t(%f %f %f)(%f %f %f)(%f %f %f)(%f %f %f)\n",
//			validBools[0], validBools[1], validBools[2], validBools[3], validBools[4], validBools[5], validBools[6], msg->x[0], msg->y[0],
//			msg->z[0], msg->x[1], msg->y[1], msg->z[1], msg->x[2], msg->y[2], msg->z[2], msg->x[3], msg->y[3], msg->z[3], msg->x[4],
//			msg->y[4], msg->z[4], msg->x[5], msg->y[5], msg->z[5], msg->x[6], msg->y[6], msg->z[6]);
//	fprintf(stdout, "uta::attitudeEstimator::Matrix3Xd getTestConstellation%d(){"
//			"	uta::attitudeEstimator::Matrix3Xd ret(3, 7);\n"
//			"	ret << %f, %f, %f, %f, %f, %f, %f, //\n"
//			"	%f, %f, %f, %f, %f, %f, %f, //\n"
//			"	%f, %f, %f, %f, %f, %f, %f;"
//			"	return ret;}\n", count, msg->x[0], msg->x[1], msg->x[2], msg->x[3], msg->x[4], msg->x[5], msg->x[6], msg->y[0], msg->y[1],
//			msg->y[2], msg->y[3], msg->y[4], msg->y[5], msg->y[6], msg->z[0], msg->z[2], msg->z[2], msg->z[3], msg->z[4], msg->z[5],
//			msg->z[6]);
//}
void HumePoseKalmanFilter::updateSelector(const Matrix<bool,1,7> &valid_leds)
{
	selector.setZero(7, 7);
	numberOfPoints = 0;
	for (int i = 0; i < 7; i++)
		if (valid_leds[0,i])
			selector(i, numberOfPoints++) = 1.0;
}
void HumePoseKalmanFilter::updateSVDMetaData(int32_t validBits)
{
	double averageVisibleLEDSAlphaFilterValue = 0.5;
	averageVisibleLEDs = averageVisibleLEDs * (1 - averageVisibleLEDSAlphaFilterValue)
			+ numberOfPoints * averageVisibleLEDSAlphaFilterValue;
	averageSingularValuesOfConstellation = averageSingularValuesOfConstellation * (1 - averageVisibleLEDSAlphaFilterValue)
			+ averageVisibleLEDSAlphaFilterValue * constellationSvds[validBits];
}
SpatialTransform<worldFrame, ledFrame> HumePoseKalmanFilter::solveLEDSubproblem(
		const SpatialTransform<worldFrame, ledFrame>& estimatedHistoricLEDPoseInWorldFrame,
		const Matrix<double,3,7>& led_mesaurements,const Matrix<bool,1,7>& led_selector)
{
	Matrix3Xd ledsInCamera(led_mesaurements);
	ledsInCamera *= 1.0 / 1000.0;
//	if (count++ % 1000 == 0)
//		printMessage();
	SpatialTransform<worldFrame, ledFrame> newEstimatedHistoricLEDPoseInWorldFrame;
	updateSelector(led_selector);
	int32_t validBits{0};
	for (int i=0;i<7;i++)
		validBits+= led_selector(0,i)*1<<i;
	if (alternateEstimators.count(validBits) == 0)
	{
		gu::attitudeEstimator::Matrix3Xd selectedConstellation = ledsInLEDFrame * selector.block(0, 0, 7, numberOfPoints);
		constellationSvds[validBits] = gu::attitudeEstimator::svdConstellation(selectedConstellation);
		alternateEstimators[validBits] = gu::attitudeEstimator::AttitudeEstimator();
		alternateEstimators[validBits].setConstellationWithPrior(selectedConstellation, priorPointWeight,
				priorMatrixWeight);
	}
	updateSVDMetaData(validBits);
	Matrix3Xd ledsinWorld = cameraPose * ledsInCamera;
	Matrix3Xd selectedData = ledsinWorld * selector.block(0, 0, 7, numberOfPoints);
	AttitudeEstimator& estimatorToUse = alternateEstimators[validBits];

	estimatorToUse.solveDirectlyWithPrior(selectedData, estimatedHistoricLEDPoseInWorldFrame.getPosition(),
			estimatedHistoricLEDPoseInWorldFrame.getDCM());
	newEstimatedHistoricLEDPoseInWorldFrame.setPosition(estimatorToUse.getPoint());
	newEstimatedHistoricLEDPoseInWorldFrame.setOrientation(estimatorToUse.getQuaternion());

	return newEstimatedHistoricLEDPoseInWorldFrame;
}
char newBodyFrameName[] = "New Body Frame";
typedef Frame<newBodyFrameName> newBodyFrame;

void HumePoseKalmanFilter::integrateIMU(const IMUDataPoint& point)
{
	SpatialTransform<worldFrame, bodyFrame> bodyPose(this->poseHistory.back().pose);
	FrameVector<worldFrame> velocity(this->poseHistory.back().velocity);
	SpatialTransform<bodyFrame, newBodyFrame> bodyToBodyUpdate;
	FrameVector<bodyFrame> bodyVelocityInBody = bodyPose.inverse() * velocity;
	FrameVector<imuFrame> isoUnitsIMUAccel = point.acceleration.scale(9.81);
	FrameVector<bodyFrame> bodyAccelerationInBody = imuPose * isoUnitsIMUAccel + bodyPose.inverse() * inertialFrameAcceleration;
	Vector3d vector = timestep * bodyVelocityInBody.getVector() ;//+ 0.5 * timestep * timestep * bodyAccelerationInBody.getVector();
	bodyToBodyUpdate.setPosition(vector);

	FrameVector<worldFrame> velocityAfterAcceleration(velocity + bodyPose * bodyAccelerationInBody.scale(timestep));

	FrameVector<bodyFrame> angularVelocityInBody = imuPose * point.angularVelocity;
	double angleOfRotation = angularVelocityInBody.getVector().norm() * timestep;
	Vector3d axisOfRotation = angularVelocityInBody.getVector().normalized();
	if (std::isnan(axisOfRotation.norm()))
		axisOfRotation << 1.0, 0.0, 0.0; // handle case of negligable rotation

	bodyToBodyUpdate.setOrientation(versor(axisOfRotation, angleOfRotation));
	SpatialTransform<worldFrame, newBodyFrame> newBodyPose (bodyPose * bodyToBodyUpdate);
	PoseHistoryElement newElement;
	newElement.pose = newBodyPose.cast<worldFrame,bodyFrame>();
	newElement.velocity = velocityAfterAcceleration;
	poseHistory.push_back(newElement);
}
void HumePoseKalmanFilter::updateIMU(double accelX, double accelY, double accelZ, double omegaX, double omegaY, double omegaZ)
{
	FrameVector<imuFrame> accel(accelX, accelY, accelZ);
	FrameVector<imuFrame> omega(omegaX, omegaY, omegaZ);
	IMUDataPoint newPoint(accel, omega);
	imuData.push_back(newPoint); // this is thread safe
	//PUNT: auto-identify timestep.
	if (ledUpdating)
		return;	//PUNT: wait for ledUpdate to finish
	integrateIMU(newPoint);
	//PUNT: pose history should mark the pose of the imu frame, so a transform will be necessary here.
	packYawPitchRoll(yaw, pitch, roll, poseHistory.back().pose.getDCM());
	Vector3d xyz(poseHistory.back().pose.getPosition());
	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
}

}
}
