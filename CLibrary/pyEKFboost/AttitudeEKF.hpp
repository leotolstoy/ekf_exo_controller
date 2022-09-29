#ifndef ATTITUDE_EKF
#define ATTITUDE_EKF
/*
AttitudeEKF.hpp, a C++ port of AttitudeEKF.py 
 - Gray C. Thomas, Ph.D. 10/07/2020
*/
#include "AttitudeEstimatorEKFFuncs.hpp"
namespace af{
class AttitudeEKF{
    Vector6f x0;
    Matrix6f P0;
    Matrix6f F0;
    Matrix6f H;
    Matrix6f Q;
    Matrix6f R;
    Matrix3f R_gyroOnly;
    Vector6f x_state_estimate;
    Matrix6f P_covar_estimate;
    Matrix6f F;
    bool isUpdateTime;
    bool first;

    float accelNormCutoff;

public:
    AttitudeEKF();
    bool isUpdateR;
    Vector6f z_model;
    Vector6f z_measured;
    Vector6f z_residual;
    Vector6f y_residual;
    void step(int, float, bool isUpdateTime=true);
    void measure(int, Vector3f, Vector3f, bool);
    Vector3f get_euler_angles();
    std::tuple<float, float> get_useful_angles(float, int sideMultiplier = -1);
};
}//af
#endif // ATTITUDE_EKF