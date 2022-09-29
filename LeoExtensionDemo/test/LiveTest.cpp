
#include <Eigen/Dense>
#include "TestAttitudeEstimator.hpp"
#include "LowPassPositionFilter.hpp"
#include "Daemon.hpp"
#include "AggressiveTimingDaemon.hpp"
#include "MatrixReceiver.hpp"
using namespace Eigen;
using gu::attitudeEstimator::getNewTestConstellation1;
class TestingDaemon1: public PhaseSpaceUDPReceiver
{
	int count = 0;
	gu::attitudeEstimator::AttitudeEstimator estimator;
	Vector3d originShift;
protected:
	void handleUpdate()
	{
		Int32BitField* bitfield = (Int32BitField*) &mostRecentMessage.validBits;
		message* msg = &mostRecentMessage;
		count++;
		if (count % 100 == 0)
		{
			fprintf(stdout,
					"Received Data: {%d,%d,%d, %d,%d,%d,%d} \n\t(%f %f %f)(%f %f %f)(%f %f %f)\n\t(%f %f %f)(%f %f %f)(%f %f %f)(%f %f %f)\n",
					bitfield->b0, bitfield->b1, bitfield->b2, bitfield->b3, bitfield->b4, bitfield->b5, bitfield->b6,
					msg->x[0], msg->y[0], msg->z[0], msg->x[1], msg->y[1], msg->z[1], msg->x[2], msg->y[2],
					msg->z[2], msg->x[3], msg->y[3], msg->z[3], msg->x[4], msg->y[4], msg->z[4], msg->x[5],
					msg->y[5], msg->z[5], msg->x[6], msg->y[6], msg->z[6]);

			if (bitfield->b0 && bitfield->b1 && bitfield->b2 && bitfield->b3 && bitfield->b4 && bitfield->b5
					&& bitfield->b6)
			{
				printf("data is good\n");
				gu::attitudeEstimator::Matrix3Xd data(3, 7);
				for (int i = 0; i < 7; i++)
					data.col(i) << msg->x[i] - originShift(0), msg->y[i] - originShift(1), msg->z[i] - originShift(2);
				std::cout << "raw data input\n" << data << std::endl;
				estimator.solveDirectly(data);
			}
			else
				printf("data has omissions\n");
		}
	}
public:
	TestingDaemon1()
	{
		gu::attitudeEstimator::Matrix3Xd constellation = getNewTestConstellation1();
		originShift = constellation.rowwise().sum() * (1.0 / 7.0);
		constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);
		estimator.setConstellation(constellation);
	}
};
class TestingDaemon2: public PhaseSpaceUDPReceiver
{
	int count = 0;
	gu::attitudeEstimator::AttitudeEstimator defaultEstimator;
	std::map<uint32_t, gu::attitudeEstimator::AttitudeEstimator> alternateEstimators;
	Vector3d originShift;
	Matrix<double, 7, 7> selector;
	unsigned numberOfPoints = 0;
protected:
	void lazyInitializeAttitudeEstimator()
	{
		if (alternateEstimators.count(mostRecentMessage.validBits) == 0)
		{
			gu::attitudeEstimator::Matrix3Xd constellation = getNewTestConstellation1();
			originShift = constellation.rowwise().sum() * (1.0 / 7.0);
			constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);
			gu::attitudeEstimator::Matrix3Xd selectedConstellation = constellation
					* selector.block(0, 0, 7, numberOfPoints);
			alternateEstimators[mostRecentMessage.validBits] = gu::attitudeEstimator::AttitudeEstimator();
			alternateEstimators[mostRecentMessage.validBits].setConstellation(selectedConstellation);
		}
	}
	void handleUpdate()
	{
		Int32BitField* bitfield = (Int32BitField*) &mostRecentMessage.validBits;
		bool validBools[7] =
		{ (bool) bitfield->b0, (bool) bitfield->b1, (bool) bitfield->b2, (bool) bitfield->b3,
			(bool) bitfield->b4, (bool) bitfield->b5, (bool) bitfield->b6 };
		message* msg = &mostRecentMessage;
		count++;
		if (count % 200 == 0)
		{
			fprintf(stdout,
					"Received Data: {%d,%d,%d, %d,%d,%d,%d} \n\t(%f %f %f)(%f %f %f)(%f %f %f)\n\t(%f %f %f)(%f %f %f)(%f %f %f)(%f %f %f)\n",
					validBools[0], validBools[1], validBools[2], validBools[3], validBools[4], validBools[5],
					validBools[6], msg->x[0], msg->y[0], msg->z[0], msg->x[1], msg->y[1], msg->z[1], msg->x[2],
					msg->y[2], msg->z[2], msg->x[3], msg->y[3], msg->z[3], msg->x[4], msg->y[4], msg->z[4],
					msg->x[5], msg->y[5], msg->z[5], msg->x[6], msg->y[6], msg->z[6]);
			fprintf(stdout, "gu::attitudeEstimator::Matrix3Xd getTestConstellation%d(){"
					"	gu::attitudeEstimator::Matrix3Xd ret(3, 7);\n"
					"	ret << %f, %f, %f, %f, %f, %f, %f, //\n"
					"	%f, %f, %f, %f, %f, %f, %f, //\n"
					"	%f, %f, %f, %f, %f, %f, %f;"
					"	return ret;}\n", count, msg->x[0], msg->x[1], msg->x[2], msg->x[3], msg->x[4], msg->x[5],
					msg->x[6], msg->y[0], msg->y[1], msg->y[2], msg->y[3], msg->y[4], msg->y[5], msg->y[6],
					msg->z[0], msg->z[2], msg->z[2], msg->z[3], msg->z[4], msg->z[5], msg->z[6]);
		}
		if (bitfield->b0 && bitfield->b1 && bitfield->b2 && bitfield->b3 && bitfield->b4 && bitfield->b5
				&& bitfield->b6)
		{
			gu::attitudeEstimator::Matrix3Xd data(3, 7);
			for (int i = 0; i < 7; i++)
				data.col(i) << msg->x[i] - originShift(0), msg->y[i] - originShift(1), msg->z[i] - originShift(2);
//			std::cout << "raw data input\n" << data << std::endl;
			defaultEstimator.solveDirectly(data);
		}
		else
		{
			selector.setZero(7, 7);
			numberOfPoints = 0;
			for (int i = 0; i < 7; i++)
				if (validBools[i])
					selector(i, numberOfPoints++) = 1.0;
			if (numberOfPoints < 4)
			{
				printf("we are not ready for this few points\n");
				return;
			}
			lazyInitializeAttitudeEstimator();

			gu::attitudeEstimator::Matrix3Xd data(3, 7);
			for (int i = 0; i < 7; i++)
				data.col(i) << msg->x[i] - originShift(0), msg->y[i] - originShift(1), msg->z[i] - originShift(2);
			gu::attitudeEstimator::Matrix3Xd selectedData = data * selector.block(0, 0, 7, numberOfPoints);
//			std::cout << "selected input data \n" << selectedData << std::endl;
			alternateEstimators[mostRecentMessage.validBits].solveDirectly(selectedData);
		}
	}
public:
	TestingDaemon2()
	{
		gu::attitudeEstimator::Matrix3Xd constellation = getNewTestConstellation1();
		originShift = constellation.rowwise().sum() * (1.0 / 7.0);
		constellation -= originShift.asDiagonal() * gu::attitudeEstimator::Matrix3Xd::Ones(3, 7);
		defaultEstimator.setConstellation(constellation);
	}
};

class AttitudeYellingDaemon: public gu::timing::AggressiveTimingDaemon
{
	gu::attitudeEstimator::LowPassPositionFilter& filter;

public:
	AttitudeYellingDaemon(gu::attitudeEstimator::LowPassPositionFilter& filter, double period) :
				AggressiveTimingDaemon(period),
				filter(filter)
	{
	}
protected:
	void agressivelyTimedLoop()
	{
		Vector3d vec = filter.getPosition();
		Quaterniond quat = filter.getOrientation();
		double leds = filter.getAverageNumberOfLEDs();
		Vector3d svds = filter.getAverageSVDForOrientationEstimation();

		fprintf(stdout,
				"timing error = %3d @ %1.3f LEDs, % 8.3E, % 8.3E, % 8.3E = [% 7.2f, % 7.2f, % 7.2f -- % 7.4f, % 7.4f, % 7.4f, % 7.4f]\n",
				getMicroError(), leds, svds(0), svds(1), svds(2), vec(0), vec(1), vec(2), quat.w(), quat.x(),
				quat.y(), quat.z());

	}
};
void testLive(void)
{
	gu::attitudeEstimator::LowPassPositionFilter tester(7 * 25, 7.881E+04 * 25);
	gu::humeEstimation::MatrixReceiver<7> receiver;
	receiver.add_listener(&tester);
	AttitudeYellingDaemon yeller(tester, 0.5);
	receiver.start();
	yeller.start();
	pause();
	puts("Terminated.");
}
