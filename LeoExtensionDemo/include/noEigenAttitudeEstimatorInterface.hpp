#ifndef NO_EIGEN_ATTITUDE_ESTIMATOR_INTERFACE
#define NO_EIGEN_ATTITUDE_ESTIMATOR_INTERFACE
namespace gu
{
namespace attitudeEstimator
{
namespace noEigenInterface
{
class HumePoseEstimator
{
public:
	virtual ~HumePoseEstimator()
	{
	}
	virtual double getPositionX()=0;
	virtual double getPositionY()=0;
	virtual double getPositionZ()=0;
	virtual double getRoll()=0;
	virtual double getYaw()=0;
	virtual double getPitch()=0;
	virtual void updateIMU(double accelX, double accelY, double accelZ, double omegaX, double omegaY,
							double omegaZ)=0;
};
HumePoseEstimator* getDefaultPoseEstimator();
HumePoseEstimator* getLowPassEstimator();
}
}
}
#endif
