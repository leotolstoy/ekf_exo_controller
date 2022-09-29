#ifndef LOW_PASS_POSITION_FILTER__HPP
#define LOW_PASS_POSITION_FILTER__HPP

#include "Eigen/Dense"
#include "MatrixListener.hpp"
#include "MatrixReceiver.hpp"
#include "AttitudeEstimator.hpp"
#include "noEigenAttitudeEstimatorInterface.hpp"
#include <map>

using namespace Eigen;
namespace gu
{
namespace attitudeEstimator
{
class LowPassPositionFilter: public gu::humeEstimation::MatrixListener<7>,
		public noEigenInterface::HumePoseEstimator
{
	int count = 0;
	gu::attitudeEstimator::AttitudeEstimator defaultEstimator;
	std::map<uint32_t, gu::attitudeEstimator::AttitudeEstimator> alternateEstimators;
	Vector3d defaultSvd;
	std::map<uint32_t, Vector3d> constellationSvds;
	Vector3d originShift;
	Matrix<double, 7, 7> selector;
	Vector3d currentPosition;
	Matrix3d currentDirectionCosineMatrix;
	Quaterniond currentQuaternion;
	double priorPointWeight = 7.0, priorMatrixWeight = 80000.0;
	unsigned numberOfPoints = 0;
	double averageVisibleLEDs = 7;
	Vector3d averageSingularValuesOfConstellation;
	double x, y, z, roll, pitch, yaw;
	std::vector<gu::humeEstimation::MatrixReceiver<7>*> managedReceivers;
public:
	LowPassPositionFilter(double priorPointWeight, double priorMatrixWeight);
	~LowPassPositionFilter();
	inline void add_dependent_receiver(gu::humeEstimation::MatrixReceiver<7>* receiver){managedReceivers.push_back(receiver);}
protected:
	void lazyInitializeAttitudeEstimator(int32_t validBits);
	void handleUpdate(const Matrix<double,3,7>&,const Matrix<bool,1,7>&);
public:
	const Vector3d& getPosition() const
	{
		return currentPosition;
	}
	const Quaterniond& getOrientation() const
	{
		return currentQuaternion;
	}
	const Matrix3d& getDCM() const
	{
		return currentDirectionCosineMatrix;
	}
	double getAverageNumberOfLEDs() const
	{
		return averageVisibleLEDs;
	}
	Vector3d getAverageSVDForOrientationEstimation() const
	{
		return Vector3d(averageSingularValuesOfConstellation);
	}
	virtual double getPositionX()
	{
		return x;
	}
	virtual double getPositionY()
	{
		return y;
	}
	virtual double getPositionZ()
	{
		return z;
	}
	virtual double getRoll()
	{
		return roll;
	}
	virtual double getYaw()
	{
		return yaw;
	}
	virtual double getPitch()
	{
		return pitch;
	}
	virtual void updateIMU(double accelX, double accelY, double accelZ, double omegaX, double omegaY, double omegaZ)
	{
	}
	;
};
}
}
#endif
