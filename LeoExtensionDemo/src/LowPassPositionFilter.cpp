#include "LowPassPositionFilter.hpp"
#include "AttitudeEstimator.hpp"
#include <SpatialTransform.hpp>
#include <Eigen/Dense>
#include <iostream>
namespace gu
{
namespace attitudeEstimator
{
gu::attitudeEstimator::Matrix3Xd getDefaultConstellation()
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
LowPassPositionFilter::LowPassPositionFilter(double priorPointWeight, double priorMatrixWeight) :
			defaultEstimator(),
			alternateEstimators(),
			originShift(),
			selector(),
			currentPosition(),
			currentDirectionCosineMatrix(),
			currentQuaternion(),
			x(0),
			y(0),
			z(0),
			yaw(0),
			pitch(0),
			roll(0),
			managedReceivers()
{
	this->priorPointWeight = priorPointWeight;
	this->priorMatrixWeight = priorMatrixWeight;
	gu::attitudeEstimator::Matrix3Xd constellation = getDefaultConstellation();
	originShift = constellation.rowwise().sum() * (1.0 / 7.0);
	constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);
	defaultEstimator.setConstellationWithPrior(constellation, priorPointWeight, priorMatrixWeight);
	defaultSvd = gu::attitudeEstimator::svdConstellation(constellation);
	std::cout<<"initializing LowPassPositionFilter"<<std::endl;
}
LowPassPositionFilter::~LowPassPositionFilter()
{
	for (auto it=managedReceivers.begin();it!=managedReceivers.end();++it)
	{
		delete (*it);
	}
}
void LowPassPositionFilter::lazyInitializeAttitudeEstimator(int32_t validBits)
{
	if (alternateEstimators.count(validBits) == 0)
	{
		gu::attitudeEstimator::Matrix3Xd constellation = getDefaultConstellation();
		originShift = constellation.rowwise().sum() * (1.0 / 7.0);
		constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);
		gu::attitudeEstimator::Matrix3Xd selectedConstellation = constellation
				* selector.block(0, 0, 7, numberOfPoints);
		constellationSvds[validBits] = gu::attitudeEstimator::svdConstellation(
				selectedConstellation);
		alternateEstimators[validBits] = gu::attitudeEstimator::AttitudeEstimator();
		alternateEstimators[validBits].setConstellationWithPrior(selectedConstellation,
				priorPointWeight, priorMatrixWeight);
	}
}
void LowPassPositionFilter::handleUpdate(const Matrix<double,3,7>& led_mesaurements,const Matrix<bool,1,7>& led_selector)
{
	int32_t validBits{0};
	for (int i=0;i<7;i++)
		validBits+= led_selector(0,i)*1<<i;

	selector.setZero(7, 7);
	numberOfPoints = 0;
	for (int i = 0; i < 7; i++)
		if (led_selector[0,i])
			selector(i, numberOfPoints++) = 1.0;

	double alpha = 0.5;
	averageVisibleLEDs = averageVisibleLEDs * (1 - alpha) + numberOfPoints * alpha;


	if (validBits==0b01111111)
	{
		gu::attitudeEstimator::Matrix3Xd data(3, 7);
		for (int i = 0; i < 7; i++)
			data.col(i) << led_mesaurements(0,i) - originShift(0), led_mesaurements(1,i) - originShift(1), led_mesaurements(2,i) - originShift(2);
//			std::cout << "raw data input\n" << data << std::endl;
		defaultEstimator.solveDirectlyWithPrior(data, currentPosition, currentDirectionCosineMatrix);
		currentDirectionCosineMatrix << defaultEstimator.getDCM();
		currentPosition << defaultEstimator.getPoint();
		averageSingularValuesOfConstellation = averageSingularValuesOfConstellation * (1 - alpha)
				+ alpha * defaultSvd;
	}
	else
	{
		lazyInitializeAttitudeEstimator(validBits);

		gu::attitudeEstimator::Matrix3Xd data(3, 7);
		for (int i = 0; i < 7; i++)
			data.col(i) << led_mesaurements(0,i) - originShift(0), led_mesaurements(1,i) - originShift(1), led_mesaurements(2,i) - originShift(2);
		gu::attitudeEstimator::Matrix3Xd selectedData = data * selector.block(0, 0, 7, numberOfPoints);
//			std::cout << "selected input data \n" << selectedData << std::endl;
		alternateEstimators[validBits].solveDirectlyWithPrior(selectedData, currentPosition,
				currentDirectionCosineMatrix);
		currentDirectionCosineMatrix << alternateEstimators[validBits].getDCM();
		currentPosition << alternateEstimators[validBits].getPoint();
		currentQuaternion = alternateEstimators[validBits].getQuaternion();

		averageSingularValuesOfConstellation = averageSingularValuesOfConstellation * (1 - alpha)
				+ alpha * constellationSvds[validBits];
	}

	count++;
//	if (count % 200 == 201)
//	{
//		fprintf(stdout,
//				"Received Data: {%d,%d,%d, %d,%d,%d,%d} \n\t(%f %f %f)(%f %f %f)(%f %f %f)\n\t(%f %f %f)(%f %f %f)(%f %f %f)(%f %f %f)\n",
//				validBools[0], validBools[1], validBools[2], validBools[3], validBools[4], validBools[5],
//				validBools[6], msg->x[0], msg->y[0], msg->z[0], msg->x[1], msg->y[1], msg->z[1], msg->x[2],
//				msg->y[2], msg->z[2], msg->x[3], msg->y[3], msg->z[3], msg->x[4], msg->y[4], msg->z[4], msg->x[5],
//				msg->y[5], msg->z[5], msg->x[6], msg->y[6], msg->z[6]);
//		fprintf(stdout, "uta::attitudeEstimator::Matrix3Xd getTestConstellation%d(){"
//				"	uta::attitudeEstimator::Matrix3Xd ret(3, 7);\n"
//				"	ret << %f, %f, %f, %f, %f, %f, %f, //\n"
//				"	%f, %f, %f, %f, %f, %f, %f, //\n"
//				"	%f, %f, %f, %f, %f, %f, %f;"
//				"	return ret;}\n", count, msg->x[0], msg->x[1], msg->x[2], msg->x[3], msg->x[4], msg->x[5],
//				msg->x[6], msg->y[0], msg->y[1], msg->y[2], msg->y[3], msg->y[4], msg->y[5], msg->y[6], msg->z[0],
//				msg->z[2], msg->z[2], msg->z[3], msg->z[4], msg->z[5], msg->z[6]);
//	}
	// noEigenInterfaceStuff
//	std::cout<<"DCM\n"<<getDCM()<<std::endl;
	Matrix3d cameraDCM;
	cameraDCM << 0, 0, -1, -1, 0, 0, 0, 1, 0;
	Vector3d currentPositionInWorld=cameraDCM*getPosition();
	x = currentPositionInWorld[0];
	y = currentPositionInWorld[1];
	z = currentPositionInWorld[2];
//	std::cout<<"cameraDCM*getDCM()*cameraDCM.inverse()\n"<<cameraDCM*getDCM()*cameraDCM.inverse()<<std::endl;
	Matrix3d res=cameraDCM*getDCM()*cameraDCM.inverse();
	gu::vectorGeometry::packYawPitchRoll(yaw, pitch, roll, res);
}

namespace noEigenInterface
{
HumePoseEstimator* getLowPassEstimator()
{
	LowPassPositionFilter* filter = new LowPassPositionFilter(7,80000);
	gu::humeEstimation::MatrixReceiver<7> *receiver= new gu::humeEstimation::MatrixReceiver<7>();
	receiver->add_listener(filter);
	filter->add_dependent_receiver(receiver);// Set proper destructor for the receiver.
	receiver->start();
	return filter;
}
}
}
}
