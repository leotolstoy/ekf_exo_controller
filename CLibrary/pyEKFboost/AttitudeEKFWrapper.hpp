#ifndef ATTITUDEEKFWRAPPER_HPP
#define ATTITUDEEKFWRAPPER_HPP

#include "PyExt.hpp"
#include "AttitudeEstimatorEKFFuncs.hpp"
#include "AttitudeEKF.hpp"

namespace gupy = gu::py_ext;
namespace py = boost::python;

class AttitudeEKFWrapper: public af:AttitudeEKF{
    // Vector6f x0;
    // Matrix6f P0;
    // Matrix6f F0;
    // Matrix6f H;
    // Matrix6f Q;
    // Matrix6f R;
    // Matrix3f R_gyroOnly;
    // Vector6f x_state_estimate;
    // Matrix6f P_covar_estimate;
    // Matrix6f F;
    // bool isUpdateTime;
    // bool first;

    // float accelNormCutoff;

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





#endif